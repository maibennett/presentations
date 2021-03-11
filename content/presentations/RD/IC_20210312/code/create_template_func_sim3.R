
create_template <- function(data, id, running_var, cov1, cov2, bandwidth, size = 1000, 
                            nsamples = 10, seed = 1000){
  
  set.seed(seed)
  
  rv = which(names(data)==running_var)
  
  data = as.data.frame(data)
  
  d = data[(data[,rv]>= bandwidth[1]) & (data[,rv] <= bandwidth[2]),]
  d_template = d
  
  sample<-rep(NA,ncol(d)+1)
  
  #number of samples
  n_samples <- nsamples
  #Sample size
  n_size <- size
  
  for(i in 1:n_samples){
    d_aux<-d
    aux <- d_aux[sample(1:nrow(d_aux), n_size, replace=FALSE),]
    aux<-cbind(aux,i)
    sample<-rbind(sample,aux)
  }
  
  sample<-sample[-1,]
  
  #Now we assess which sample is more representative:
  
  #Most important covariates are: cov1
  #Other cov: cov2
  
  d_covariates<-as.data.frame(d[,c(cov1,cov2)])

  colnames(d_covariates)<- c(cov1,cov2)
    
  d_mean<-apply(d_covariates,MARGIN=2,FUN=mean,na.rm=1)
  
  S<-var(d_covariates,na.rm=TRUE)
  
  s_covariates<-as.data.frame(cbind(sample[,c(cov1,cov2)],sample$i))
  
  colnames(s_covariates)<-c(cov1,cov2,"group")
  
  s_mean<-matrix(NA,nrow=n_samples,ncol=ncol(s_covariates))
  
  for (j in 1:n_samples){
    s_mean[j,]<-apply(s_covariates[s_covariates$group==j,],MARGIN=2,FUN=mean,na.rm=1)
  }
  
  
  #Mahalanobis distance
  
  #We create an empty vector to store the sum of the mahalanobis distance for the covariates
  dist<-rep(NA,n_samples)
  min<-100000
  min_index<-0
  
  #to avoid singularity problem
  suppressMessages(library(Smisc))
  lines = findDepMat(S, rows = FALSE, tol = 1e-17)
  colnames(S)[which(lines==TRUE)]
  
  S_inv = solve(S, tol = 1e-28)
  
  for (i in 1:n_samples){
    #Incorporate adjustment for redundancy of highly correlated variables
    #dist[i]<-mahalanobis(s_mean[i,1:(ncol(s_mean)-1)],d_mean,S_inv, inverted = TRUE)
    dist[i]<-mahalanobis(s_mean[i,1:(ncol(s_mean)-1)],d_mean,S)
    
    #We update the minimum of diff (and its index) in every loop
    if(min(dist[!is.na(dist)])<min){
      min<-min(dist[!is.na(dist)])
      min_index<-i
    }
    #print(i)
  }
  
  #The most representative of the samples obtained is:
  min_index
  
  #Compare the sample with the complete population:
  sample <- as.data.frame(sample)
  
  template_id<-sample[sample$i==min_index,id]
  template<-d[d[,id] %in% template_id,]
  
  comparison<-cbind(round(d_mean,2),round(s_mean[min_index,1:(ncol(s_mean)-1)],2))
  
  suppressMessages(library("xtable"))
  
  rownames(comparison)<-c(cov1,cov2)
  
  #comparison_latex<-xtable(comparison)
  
  #print(comparison_latex,floating=F,include.rownames=T,include.colnames=T,file="template_comparison_sim.tex")
  
  return(list(template = template, d_template = d_template, comparison = comparison, bandwidth = bandwidth))
}