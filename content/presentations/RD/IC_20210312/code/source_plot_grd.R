## Self contained code

# Title: Simulations
# Paper: How far is too far?
# Authors: M. Bennett
# Date: 08-06-19

## Load libraries
library(designmatch)
library(Rglpk)
library(gurobi)
library(xtable)
library(dplyr)
library(readstata13)
library(rdrobust)
library(exact2x2)
library(simstudy)
library(exact2x2)
library(foreign)
library(rbounds)
library(stats)
library(sensitivitymv)
library(DTComPair)
library(Rlab)
library(nprobust)
library(KernSmooth)
library(lpdensity)
###################################################
#### Load real data ###############################
###################################################

#Type of simulation

#Correlation: high vs low (covariate-running variable)
correlation = "low"

#Sample size: small vs large
sample_size = "large"

#relationship between outcome and running variable: (linear vs cubic)
rel_y_r = "linear"

#Unobserved confounder: (none vs linear)
unobs = "linear"

#####################################################################
######## Parameters #################################################
#####################################################################

# Generalization bandwidth:
b = 200

# N simulations
sim_start = 1
sim_end = 500

seeds = seq(sim_start,sim_end,1)

results = rep(NA,20)
lower.ci = rep(NA,8)
upper.ci = rep(NA,8)

sim = 33

print("###############################################")
print("###############################################")
print("###############################################")
print(paste("########### SIM:",sim,"#################################"))
print("###############################################")
print("###############################################")
set.seed(sim)

if(sample_size=="large"){n = 20000}
if(sample_size=="small"){n = 2000}

# Define observable covariate:
x0 = rnorm(n,0,10)
# Define unobserved covariate:
u0 = rnorm(n,0,10)

# Create the outcome:
if(rel_y_r=="linear"){
  # Define running variable:
  r0 = 5*x0 + 10*u0 + rnorm(n,0,100)
  
  r0[r0>-b & r0<b] = 5*x0[r0>-b & r0<b] + rnorm(sum(r0>-b & r0<b),0,140) 
  
  #y0 = 0.02*r0 + rnorm(n,0,10)
  
  y0 = x0 + u0 + rnorm(n,0,10)
  y0[r0>-b & r0<b] = 2*x0[r0>-b & r0<b] + rnorm(sum(r0>-b & r0<b),0,10)
}
if(rel_y_r=="cubic"){
  # Define running variable:
  r0 = 5*x0 + 10*u0 + rnorm(n,0,100)
  
  r0[r0>-b & r0<b] = 1*x0[r0>-b & r0<b] - 0.01*(x0[r0>-b & r0<b])^2 + 0.001*(x0[r0>-b & r0<b])^3 + 
    rnorm(sum(r0>-b & r0<b),0,140) 
  
  y0 = 0.002*r0 - 0.00001*r0^2 + 0.0000001*r0^3 + rnorm(n,0,10)
  
  #y0 = 1*x0 - 0.01*x0^2 + 0.001*x0^3 + u0 + rnorm(n,0,10)
  #y0[r0>-b & r0<b] = 4*x0[r0>-b & r0<b] + - 0.1*(x0[r0>-b & r0<b])^2 + 0.1*(x0[r0>-b & r0<b])^3 + 
  #  rnorm(sum(r0>-b & r0<b),0,10)
}


# Create the outcome:
if(correlation=="high"){
  # Define running variable:
  r0 = 10*x0 + 5*u0 + rnorm(n,0,100)
  
  r0[r0>-b & r0<b] = 10*x0[r0>-b & r0<b] + rnorm(sum(r0>-b & r0<b),0,140) 
  
  #y0 = 0.02*r0 + rnorm(n,0,10)
  
  y0 = x0 + u0 + rnorm(n,0,10)
  y0[r0>-b & r0<b] = 2*x0[r0>-b & r0<b] + rnorm(sum(r0>-b & r0<b),0,10)
}


dt0 = as.data.frame(cbind(r0,x0,u0,y0))
names(dt0) = c("r","x","u","y")

threshold0 = mean(dt0$r)

dt0$running_var = dt0$r - threshold0

rdplot(dt0$y,dt0$running_var,title="RD Plot for T=0",x.label = "Running variable",
       y.label = "Outcome")

#Post-intervention data

# Define observable covariate:
x1 = rnorm(n,0,10)
# Define unobserved covariate:
u1 = rnorm(n,0,10)

r1 = 5*x1 + 10*u1 + rnorm(n,0,100)

# Create the outcome:
if(rel_y_r=="linear"){
  r1[r1>-b & r1<b] = 5*x1[r1>-b & r1<b] + rnorm(sum(r1>-b & r1<b),0,140)
  
  #y1 = 0.02*r1 + rnorm(n,0,10)
  
  y1 = x1 + u1
  y1[r1>-b & r1<b] = 2*x1[r1>-b & r1<b]
}

# Create the outcome:
if(correlation=="high"){
  # Define running variable:
  r1 = 10*x1 + 5*u1 + rnorm(n,0,100)
  
  r1[r1>-b & r1<b] = 10*x1[r1>-b & r1<b] + rnorm(sum(r1>-b & r1<b),0,140) 
  
  y1 = x1 + u1
  y1[r1>-b & r1<b] = 2*x1[r1>-b & r1<b]
}

dt1 = as.data.frame(cbind(r1,x1,u1,y1))
names(dt1) = c("r","x","u","y")

threshold1 = mean(dt1$r)

dt1$running_var = dt1$r - threshold1

dt0$t = 0
dt1$t = 1

dt0$treat=0
dt0$treat[dt0$r< threshold0] = 1

dt1$treat=0
dt1$treat[dt1$r< threshold1] = 1

#Transform important covariates to quintiles:
cov1 = c("x")
cov2 = c()

dt0 = as.data.frame(dt0)
dt1 = as.data.frame(dt1)

cov1_10_0 = list()
cov1_10_1 = list()

for(l in 1:length(cov1)){
  
  names0 = names(dt0)
  
  q0 = quantile(dt0[,cov1[l]], probs = seq(0, 1, by = .1))
  
  cov1_10_0[[l]] = rep(10,nrow(dt0))
  cov1_10_0[[l]][dt0[,cov1[l]]<=q0[2]] = 1 
  cov1_10_0[[l]][dt0[,cov1[l]]>q0[2] & dt0[,cov1[l]]<=q0[3]] = 2
  cov1_10_0[[l]][dt0[,cov1[l]]>q0[3] & dt0[,cov1[l]]<=q0[4]] = 3
  cov1_10_0[[l]][dt0[,cov1[l]]>q0[4] & dt0[,cov1[l]]<=q0[5]] = 4
  cov1_10_0[[l]][dt0[,cov1[l]]>q0[5] & dt0[,cov1[l]]<=q0[6]] = 5
  cov1_10_0[[l]][dt0[,cov1[l]]>q0[6] & dt0[,cov1[l]]<=q0[7]] = 6
  cov1_10_0[[l]][dt0[,cov1[l]]>q0[7] & dt0[,cov1[l]]<=q0[8]] = 7
  cov1_10_0[[l]][dt0[,cov1[l]]>q0[8] & dt0[,cov1[l]]<=q0[9]] = 8
  cov1_10_0[[l]][dt0[,cov1[l]]>q0[9] & dt0[,cov1[l]]<=q0[10]] = 9
  
  d_cov1_10_0 <- model.matrix(~ as.factor(cov1_10_0[[l]]))
  d_cov1_10_0 <- d_cov1_10_0[,-1]
  
  dt0 <- cbind(dt0,d_cov1_10_0)
  names(dt0) <- c(names0,paste0(cov1[l],"_",seq(2,10)))
  
  names1 <- names(dt1)
  
  q1 = quantile(dt1[,cov1[l]], probs = seq(0, 1, by = .1))
  
  cov1_10_1[[l]] = rep(10,nrow(dt1))
  cov1_10_1[[l]][dt1[,cov1[l]]<q1[2]] = 1 
  cov1_10_1[[l]][dt1[,cov1[l]]>q1[2] & dt1[,cov1[l]]<=q1[3]] = 2
  cov1_10_1[[l]][dt1[,cov1[l]]>q1[3] & dt1[,cov1[l]]<=q1[4]] = 3
  cov1_10_1[[l]][dt1[,cov1[l]]>q1[4] & dt1[,cov1[l]]<=q1[5]] = 4
  cov1_10_1[[l]][dt1[,cov1[l]]>q1[5] & dt1[,cov1[l]]<=q1[6]] = 5
  cov1_10_1[[l]][dt1[,cov1[l]]>q1[6] & dt1[,cov1[l]]<=q1[7]] = 6
  cov1_10_1[[l]][dt1[,cov1[l]]>q1[7] & dt1[,cov1[l]]<=q1[8]] = 7
  cov1_10_1[[l]][dt1[,cov1[l]]>q1[8] & dt1[,cov1[l]]<=q1[9]] = 8
  cov1_10_1[[l]][dt1[,cov1[l]]>q1[9] & dt1[,cov1[l]]<=q1[10]] = 9
  
  d_cov1_10_1 <- model.matrix(~ as.factor(cov1_10_1[[l]]))
  d_cov1_10_1 <- d_cov1_10_1[,-1]
  
  dt1 <- cbind(dt1,d_cov1_10_1)
  names(dt1) <- c(names1,paste0(cov1[l],"_",seq(2,10)))
}

dt = rbind(dt0,dt1)


# Matching

data = dt
tvar = "t"
grid = c(quantile(dt$running_var[dt$t==1 & dt$treat==1],probs=seq(0,1,0.2))[-6],0,
         quantile(dt$running_var[dt$t==1 & dt$treat==0],probs=seq(0,1,0.2))[-1])
bandwidth = c(max(grid[grid<0]),0)
cov1 = c(paste0(cov1[1],"_",seq(2,10)))
cov2 = "x"
out_var = "y"
tols = c(0.0001,0.05)
running_var = "running_var"
threshold = 0
treat = "treat"
size = 1000
if(sample_size=="small"){size = 100} 
nsamples = 50
seed = sim
dir_data = dir_data
t_max_alloc = 30
significance = 0.05
match_prop = 0.99
if(sample_size=="small"){match_prop = 0.95}
drop_tail=0
if(sample_size=="small"){drop_tail = 0}
print = FALSE
if(sim==1){print=TRUE}
name_plot = paste0("BW_corr_",correlation,"_size_",sample_size,"_relyr_",rel_y_r,"_unobs_",unobs)

############ grd_function_sim3.R

source(paste0(dir_data,"/functions/create_template_func_sim3.R"))

bandwidth_aux = bandwidth
bandwidth_opt = bandwidth

grid_t = grid[grid<0]
grid_t = grid_t[order(grid[grid<0],decreasing = TRUE)]

optimal_bandwidth = NA

STOP = FALSE
OPT = FALSE
at_least_one = FALSE
J = 1

data = as.data.frame(data)

data$id_grd = seq(1,nrow(data),1)

data$template = 0

sig_dif = NA

ct = create_template(data = data, id = "id_grd", running_var = running_var,
                     cov1 = cov1, cov2 = cov2,
                     bandwidth = bandwidth_aux, size = size, 
                     seed = seed)

template = ct$template

for(g in 1:(length(grid)-1)){
  
  
  print("###############################################")
  print("###############################################")
  print(paste("######### Grid Num:",g,"(",round(grid[g],2),"-",round(grid[g+1],2),")"))
  print("###############################################")
  print("###############################################")
  
  bw_data = c(grid[g],grid[g+1])
  
  d2 = data[(data[,running_var]>=bw_data[1]) & (data[,running_var]<bw_data[2]),]
  
  data_bef <- d2[d2[,tvar]==0,]
  
  template[,treat] = NA
  template[,tvar] = NA
  template$template = 1
  
  d_bef <- rbind(data_bef,template)
  
  #Generate id for each obs
  d_bef$id=seq(1,nrow(d_bef),1)
  
  if(g==1){
    d_match_bef = rep(NA, ncol(d_bef)+7)
    count0 <- 0
  }
  
  
  print("###############################################")
  print("###############################################")
  print(paste("######### MATCHING #########"))
  print("###############################################")
  print("###############################################")
  
  names_mom <- c(cov1, cov2)
  
  ##############################################
  # Match
  ##############################################
  
  d_aux <- d_bef
  
  dim(d_aux)
  
  d_aux$t_ind <- d_aux$template
  
  #We order data so we have treated observations first.
  d_aux = d_aux[order(-d_aux$t_ind), ]
  
  t_ind=d_aux$t_ind
  
  table(t_ind)
  
  mom_covs <- d_aux[,c(names_mom)]
  
  if(length(cov1)==1){mom_tols1 <- round(absstddif(cbind(1,mom_covs[,1]),t_ind,tols[1]),5)[2]}
  if(length(cov1)>1){mom_tols1 <- round(absstddif(mom_covs[,1:length(cov1)],t_ind,tols[1]),5)}
  
  if(length(cov2)==1){mom_tols2 <- round(absstddif(cbind(1,mom_covs[,ncol(mom_covs)]),t_ind,tols[2]),5)[2]}
  if(length(cov2)>1){mom_tols2 <- round(absstddif(mom_covs[,(length(cov1)+1):ncol(mom_covs)],t_ind,tols[2]),5)}
  
  mom_tols <- c(mom_tols1,mom_tols2)
  mom <- list(covs = mom_covs, tols = mom_tols)
  
  # Solver options
  t_max <- 60*t_max_alloc
  solver <- "gurobi"
  #solver <- "glpk"
  approximate <- 0
  solver <- list(name = solver, t_max = t_max, approximate = approximate, round_cplex = 0, trace_cplex = 0)
  
  # Match
  out = cardmatch(t_ind = t_ind, mom = mom, solver = solver) 
  
  t_id = out$t_id  
  c_id = out$c_id
  
  t_id_d=d_aux$id[t_id]
  c_id_d=d_aux$id[c_id]
  
  group_id = out$group_id+count0
  
  if(length(t_id)>length(c_id)){
    
    diff = seq(length(c_id)+1,length(t_id),1)
    t_id = out$t_id[-diff] 
    
    t_id_d=d_aux$id[t_id]
    
    group_id = out$group_id[-diff]+count0
  }
  
  if(length(t_id)<length(c_id)){
    
    diff = seq(length(t_id)+1,length(c_id),1)
    c_id = out$c_id[-diff] 
    
    c_id_d=d_aux$id[c_id]
    
    group_id = out$group_id[-diff]+count0
  }
  
  
  if (length(t_id)>0 & length(t_id)==length(c_id)) {
    cycle = FALSE
    
    prop = length(c_id)/min(table(t_ind))
    
    d_aux_2 = cbind(d_aux[c(t_id, c_id), ],t_id,c_id, 
                    group_id,
                    t_id_d,c_id_d,g,(g+1), prop)
    d_match_bef = rbind(d_match_bef, d_aux_2)
    if (count0 == 0) {
      d_match_bef = d_match_bef[-1, ]
    }
    count0 = nrow(d_match_bef)/2
  }
  
  cat("\n", "*************************************************", "\n", sep = "")
  cat("\n", "* Matching Group Number: ",1, sep = "")
  cat("\n", "* Matching Group and number of observations: ", sep = "")
  cat("\n", "* Number of matched pairs: ", length(c_id), sep = "")
  cat("\n", "* Proportion of possible pairs matched: ",round(length(c_id)/min(table(t_ind)), 3), sep="" )
  cat("\n", "* Matching time (mins): ", round(out$time/60, 2), sep = "")
  cat("\n", "*************************************************", "\n", sep = "")
  
  #### Re-matching
  
  if(g==1){
    d_rematch_bef = rep(NA, ncol(d_match_bef)+5)
    count0_r <- 0
  }
  
  
  d_aux <- d_aux_2
  
  d_aux$t_ind <- d_aux$template
  
  table(d_aux$t_ind)
  
  #We order data so we have treated observations first.
  d_aux = d_aux[order(-d_aux$t_ind), ]
  #attach(d_aux)
  #names(d_aux)
  
  t_ind=d_aux$t_ind
  
  table(t_ind)
  
  X_mat_aux <- d_aux[,c(cov1,cov2)]
  
  dist_mat <- distmat(t_ind,X_mat_aux)
  
  subset_weight <- NULL
  
  total_pairs <- sum(t_ind)
  
  n_matches <- 1
  
  # Solver options
  t_max <- 60*t_max_alloc
  solver <- "gurobi"
  approximate <- 0
  solver <- list(name = solver, t_max = t_max, approximate = approximate, round_cplex = 0, trace_cplex = 0)
  
  # Match
  out = bmatch(t_ind = t_ind, subset_weight = subset_weight, dist_mat = dist_mat, 
               total_groups = total_pairs, solver = solver) 
  
  t_idr = out$t_id  
  c_idr = out$c_id
  
  t_id_dr=d_aux$id_grd[t_idr]
  c_id_dr=d_aux$id_grd[c_idr]
  
  if (length(t_idr)>0) {
    cycle = FALSE
    group_idr = out$group_id+count0_r
    #d_aux_2 = cbind(d_aux[c(t_id, c_id), ], group_id)
    d_aux_2r = cbind(d_aux[c(t_idr, c_idr), ],t_idr,c_idr, group_idr,t_id_dr,c_id_dr)
    d_rematch_bef = rbind(d_rematch_bef, d_aux_2r)
    if (count0_r == 0) {
      d_rematch_bef = d_rematch_bef[-1, ]
    }
    count0_r = nrow(d_rematch_bef)/2
  }
}

d_match = d_rematch_bef[d_rematch_bef$template==0,]

d_match$bw_l = grid[d_match$g]
d_match$bw_u = grid[d_match$`(g + 1)`]
d_match$bw = 0.5*(d_match$bw_l + d_match$bw_u)

#Plot proportion of matched
plot(d_match$bw,d_match$prop)

#Drop matching iterations which couldn't be properly matched (missing >5%):
drop = c()

for(i in 1:length(unique(d_match$bw))){
  
  d_aux = d_match[d_match$bw==unique(d_match$bw)[i],]
  
  if(d_aux$prop[1]<match_prop){
    drop = c(drop,unique(d_match$bw)[i])
  }
}

# Pick only two values for dropping:
if(length(which(drop<threshold))>0 & length(which(drop>threshold))>0){
  matching_bw = c(drop[max(which(drop<threshold))],drop[min(which(drop>threshold))])
}

if(length(which(drop<threshold))>0 & length(which(drop>threshold))==0){
  matching_bw = c(drop[max(which(drop<threshold))],grid[length(grid)]+1)
}

if(length(which(drop<threshold))==0 & length(which(drop>threshold))>0){
  matching_bw = c(grid[1]-1,drop[min(which(drop>threshold))])
}

if(length(which(drop<threshold))==0 & length(which(drop>threshold))==0){
  matching_bw = c(grid[1]-1,grid[length(grid)]+1)
}

d_match = d_match[d_match[,"bw"]>matching_bw[1] & d_match[,"bw"]<matching_bw[2],]

#Trim ends:
ecdf_r = quantile(d_match[,running_var], probs = c(drop_tail/2,(1-drop_tail/2)))
low_r = ecdf_r[1]
up_r = ecdf_r[2]

d_match = d_match[d_match[,running_var]>=low_r & d_match[,running_var]<=up_r,]

#Where to evaluate de function:
eval = c(quantile(d_match[d_match[,running_var]<0,running_var],probs = seq(0,1,0.01))[-length(seq(0,1,0.01))],
         0,
         quantile(d_match[d_match[,running_var]>=0,running_var],probs = seq(0,1,0.01))[-1])

#IF THE MODEL CANNOT RUN:
possibleError <- tryCatch(
  lp.model <- lprobust(d_match[,out_var],d_match[,running_var],p=3,eval = eval),
  error=function(e) e
)

STOP = inherits(possibleError, "error")

if(STOP==FALSE){
  
  at_least_one = TRUE
  
  #lp.model <- lprobust(d_match$y_diff,d_match$bw)
  nprobust.plot(lp.model)
  
  est = lp.model$Estimate[,5]
  se = lp.model$Estimate[,7]
  eval = lp.model$Estimate[,1]
  
  
  cil = est - qnorm(1-significance/2)*se
  ciu = est + qnorm(1-significance/2)*se
  
  nl = which(cil[which(eval>=0)] > est[which(eval==0)])
  nl = c(nl[1],nl[-1][diff(nl)==1])
  if(is.na(nl[1])){nl = length(eval[eval>=0])}
  
  nu = which(ciu[which(eval<0)] < est[which(eval==0)])
  nu = c(nu[-length(nu)][diff(nu)==1],nu[length(nu)])
  if(is.na(nu[1])){nu = 1}
  
  gri = sort(c(unique(d_match$bw_l),max(d_match$bw_u)))
  
  opt_bw1_l = gri[1]
  opt_bw1_u = gri[length(gri)]
  
  opt_bw2_l = gri[1]
  opt_bw2_u = gri[length(gri)]
  
  bias_l = est[1]
  bias_u = est[eval==0]
  
  obw_grid1_l = gri[1]
  obw_grid1_u = gri[length(gri)]
  
  obw_grid2_l = gri[1]
  obw_grid2_u = gri[length(gri)]
  
  if(length(nu)>0){
    opt_bw1_l = eval[eval<0][nu[length(nu)]]
    bias_l = est[eval<0][nu[length(nu)]]
    obw_grid1_l = gri[min(which(gri/opt_bw1_l<1))]
  }
  
  if(length(nl)>0){
    opt_bw1_u = eval[eval>=0][nl[1]]
    obw_grid1_u = gri[max(which(gri/opt_bw1_u<1))]
  }
  
  if(opt_bw1_l>bandwidth_aux[1]){
    STOP = TRUE
  }
  
  if(opt_bw1_u<bandwidth_aux[2]){
    STOP = TRUE
  }
  
  bias = bias_l - bias_u
  
  optimal_bandwidth = c(max(c(opt_bw1_l,opt_bw2_l)),min(c(opt_bw1_u,opt_bw2_u)))
  obw_grid = c(max(c(obw_grid1_l,obw_grid2_l)),min(c(obw_grid1_u,obw_grid2_u)))
  template_n = template
  
  if(optimal_bandwidth[1]>=bandwidth_aux[1] | (min(data[,running_var])>=bandwidth_aux[1])){
    bandwidth_opt = c(optimal_bandwidth[1],threshold)
    OPT = TRUE
  }
  
  if(optimal_bandwidth[1]<bandwidth_aux[1]){
    bandwidth_opt = bandwidth_aux
    #bandwidth_aux = c(optimal_bandwidth[1],threshold)
    bandwidth_aux = c(grid_t[J+1],threshold)
  }
  
  if(at_least_one==1){
    STOP = FALSE
  }
}

##############################################################################################

#sim 3, 33 corr low sample large
library(ggplot2)
library(firasans)

col_text = "#23373B"

#palette v2
col1 = "#8507C4"
#col2 = "#28A4FF"
#col3 = "#FFDB00"

#palette v3
#col1 = "#8A1C7C"
#col1 = "#8A1C8D"
col2 = "#EE6352"
col3 = "#FFCF4F"

par(mfrow=c(1,1),mar=c(4,3.5,1.5,2))

set.seed(1)

random = sample(1:nrow(data),2000,replace=FALSE)

data2 = data[random,]


### Plot 1

ggplot(data2[data2[,tvar]==0,],aes(y=y,x=running_var, label=NA), fill="white")+geom_point(fill="white", color="slategrey", size=2, shape=21)+
  geom_vline(aes(xintercept=grid[1]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[2]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[3]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[4]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[5]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[6]),linetype='dashed',col="dark grey",lwd=1) +
  geom_vline(aes(xintercept=grid[8]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[9]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[10]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[11]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=0),linetype='dashed',lwd=1.5) +
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="dark grey"), size=1.5)+
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="black"), size=1.5)+
  theme_ipsum_fsc(grid="Y", plot_title_size = 20,
                  axis_text_size = 10,
                  axis_title_size=12,
                  axis_col = col_text)+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  scale_color_identity(name='Potential Outcomes',
                       breaks=c('dark grey','black'),
                       labels = c("Grid","Cutoff"),
                       guide = "legend")+
  #  annotate("pointrage",x=-600,y=41,xend=-600,yend=41,shape=21,size=2,col="dark grey") +
  #  annotate("text",x=13.05,y=8.4,label = "H*", size=5.5, colour=col_text)+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=20, colour=col_text,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12, colour=col_text),
        axis.title.y = element_text(size=20, colour=col_text,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12, colour=col_text),legend.position=c(0.1, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size=15, colour=col_text),
        legend.background = element_rect(fill=alpha("white",0),colour =alpha("white",0)),
        title = element_text(size=25, colour=col_text))+
  ylim(-60,60) + xlim(-600,600)+
  annotate("pointrange",x=-588,y=35,ymin=35,ymax=35,shape=21,size=1,col="slategrey",fill="white",alpha=1)+
  annotate("text",x=-550,y=35,label = "Observations", size=5, colour=col_text,hjust = 0,family=font_fsc)




## Plot 2

#RUN THE FIRST ITERATION

#d_match = grd$d_rematch

random_t = sample(1:nrow(d_match),600,replace=FALSE)
d_match2 = d_match[random_t,]

plot(data2[data2[,tvar]==0,running_var],data2[data2[,tvar]==0,out_var], pch="o", cex=0.5,col="dark grey",
     ylim=c(-60,60),ann=FALSE)
mtext(side = 2, text = "Outcome", line = 2)

points(d_match2[d_match2[,running_var]>=bandwidth[1] & bandwidth[2]>=d_match2[,running_var],out_var] ~ 
         d_match2[d_match2[,running_var]>=bandwidth[1] & bandwidth[2]>=d_match2[,running_var],running_var],col="deepskyblue3",pch="o",cex=0.6)

for(i in 1:length(grid)){
  abline(v=grid[i],lty=2,col="black",lwd=1)
}

abline(v=0,col="black",lty=2,lwd=2)

legend("topleft",legend = c("Data","Grid","Cutoff","Matched obs in H1"),ncol=1,
       lty=c(NA,2,2,NA),pch=c("o",NA,NA,"o"),lwd=c(NA,1,2,NA),
       col=c("darkgrey","black","black","deepskyblue3"),cex=0.8)


ggplot(data2[data2[,tvar]==0,],aes(y=y,x=running_var, label=NA), fill="white")+
  geom_point(fill="white", color="slategrey", size=2, shape=21)+
  geom_point(data=d_match2[d_match2[,running_var]>=bandwidth[1] & bandwidth[2]>=d_match2[,running_var],],
             aes(x=running_var,y=y),fill=col2, color='light grey', size=2, shape=21)+
  geom_vline(aes(xintercept=grid[1]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[2]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[3]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[4]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[5]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[6]),linetype='dashed',col="dark grey",lwd=1) +
  geom_vline(aes(xintercept=grid[8]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[9]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[10]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[11]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=0),linetype='dashed',lwd=1.5) +
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="dark grey"), size=1.5)+
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="black"), size=1.5)+
  theme_ipsum_fsc(grid="Y", plot_title_size = 20,
                  axis_text_size = 10,
                  axis_title_size=12,
                  axis_col = col_text)+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  scale_color_identity(name='Potential Outcomes',
                       breaks=c('dark grey','black'),
                       labels = c("Grid","Cutoff"),
                       guide = "legend")+
  #  annotate("pointrage",x=-600,y=41,xend=-600,yend=41,shape=21,size=2,col="dark grey") +
  #  annotate("text",x=13.05,y=8.4,label = "H*", size=5.5, colour=col_text)+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=20, colour=col_text,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12, colour=col_text),
        axis.title.y = element_text(size=20, colour=col_text,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12, colour=col_text),legend.position=c(0.1, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size=15, colour=col_text),
        legend.background = element_rect(fill=alpha("white",0),colour =alpha("white",0)),
        title = element_text(size=25, colour=col_text))+
  ylim(-60,60) + xlim(-600,600)+
  annotate("pointrange",x=-588,y=34,ymin=34,ymax=34,shape=21,size=1,col=col2,fill=col2,alpha=1)+
  annotate("text",x=-550,y=34,label = "Matched Obs", size=5, colour=col_text,hjust = 0,family=font_fsc)



ggplot(data2[data2[,tvar]==0,],aes(y=y,x=running_var, label=NA), fill="white")+
  geom_point(fill="white", color="slategrey", size=2, shape=21)+
  geom_point(data=d_match2[d_match2[,running_var]>=bandwidth[1] & bandwidth[2]>=d_match2[,running_var],],
             aes(x=running_var,y=y),fill=col2, color=col2, size=2, shape=21,alpha=0.2)+
  geom_point(data=d_match2[d_match2[,running_var]>=grid[4] & grid[5]>=d_match2[,running_var],],
             aes(x=running_var,y=y),fill=col2, color='light grey', size=2, shape=21)+
  geom_vline(aes(xintercept=grid[1]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[2]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[3]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[4]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[5]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[6]),linetype='dashed',col="dark grey",lwd=1) +
  geom_vline(aes(xintercept=grid[8]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[9]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[10]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[11]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=0),linetype='dashed',lwd=1.5) +
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="dark grey"), size=1.5)+
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="black"), size=1.5)+
  theme_ipsum_fsc(grid="Y", plot_title_size = 20,
                  axis_text_size = 10,
                  axis_title_size=12,
                  axis_col = col_text)+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  scale_color_identity(name='Potential Outcomes',
                       breaks=c('dark grey','black'),
                       labels = c("Grid","Cutoff"),
                       guide = "legend")+
  #  annotate("pointrage",x=-600,y=41,xend=-600,yend=41,shape=21,size=2,col="dark grey") +
  #  annotate("text",x=13.05,y=8.4,label = "H*", size=5.5, colour=col_text)+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=20, colour=col_text,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12, colour=col_text),
        axis.title.y = element_text(size=20, colour=col_text,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12, colour=col_text),legend.position=c(0.1, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size=15, colour=col_text),
        legend.background = element_rect(fill=alpha("white",0),colour =alpha("white",0)),
        title = element_text(size=25, colour=col_text))+
  ylim(-60,60) + xlim(-600,600)+
  annotate("pointrange",x=-588,y=34,ymin=34,ymax=34,shape=21,size=1,col=col2,fill=col2,alpha=0.2)+
  annotate("text",x=-550,y=34,label = "Matched Obs", size=5, colour=col_text,hjust = 0,family=font_fsc)




ggplot(data2[data2[,tvar]==0,],aes(y=y,x=running_var, label=NA), fill="white")+
  geom_point(fill="white", color="slategrey", size=2, shape=21)+
  geom_point(data=d_match2[d_match2[,running_var]>=bandwidth[1] & bandwidth[2]>=d_match2[,running_var],],
             aes(x=running_var,y=y),fill=col2, color=col2, size=2, shape=21,alpha=0.2)+
  geom_point(data=d_match2[d_match2[,running_var]>=grid[4] & grid[5]>=d_match2[,running_var],],
             aes(x=running_var,y=y),fill=col2, color=col2, size=2, shape=21,alpha=0.2)+
  geom_point(data=d_match2[d_match2[,running_var]>=grid[3] & grid[4]>=d_match2[,running_var],],
             aes(x=running_var,y=y),fill=col2, color='light grey', size=2, shape=21,alpha=1)+
  geom_vline(aes(xintercept=grid[1]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[2]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[3]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[4]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[5]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[6]),linetype='dashed',col="dark grey",lwd=1) +
  geom_vline(aes(xintercept=grid[8]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[9]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[10]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[11]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=0),linetype='dashed',lwd=1.5) +
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="dark grey"), size=1.5)+
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="black"), size=1.5)+
  theme_ipsum_fsc(grid="Y", plot_title_size = 20,
                  axis_text_size = 10,
                  axis_title_size=12,
                  axis_col = col_text)+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  scale_color_identity(name='Potential Outcomes',
                       breaks=c('dark grey','black'),
                       labels = c("Grid","Cutoff"),
                       guide = "legend")+
  #  annotate("pointrage",x=-600,y=41,xend=-600,yend=41,shape=21,size=2,col="dark grey") +
  #  annotate("text",x=13.05,y=8.4,label = "H*", size=5.5, colour=col_text)+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=20, colour=col_text,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12, colour=col_text),
        axis.title.y = element_text(size=20, colour=col_text,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12, colour=col_text),legend.position=c(0.1, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size=15, colour=col_text),
        legend.background = element_rect(fill=alpha("white",0),colour =alpha("white",0)),
        title = element_text(size=25, colour=col_text))+
  ylim(-60,60) + xlim(-600,600)+
  annotate("pointrange",x=-588,y=34,ymin=34,ymax=34,shape=21,size=1,col=col2,fill=col2,alpha=0.2)+
  annotate("text",x=-550,y=34,label = "Matched Obs", size=5, colour=col_text,hjust = 0,family=font_fsc)



ggplot(data2[data2[,tvar]==0,],aes(y=y,x=running_var, label=NA), fill="white")+
  geom_point(fill="white", color="slategrey", size=2, shape=21)+
  geom_point(data=d_match2[d_match2[,running_var]>=bandwidth[1] & bandwidth[2]>=d_match2[,running_var],],
             aes(x=running_var,y=y),fill=col2, color=col2, size=2, shape=21,alpha=0.2)+
  geom_point(data=d_match2[d_match2[,running_var]>=grid[4] & grid[5]>=d_match2[,running_var],],
             aes(x=running_var,y=y),fill=col2, color=col2, size=2, shape=21,alpha=0.2)+
  geom_point(data=d_match2[d_match2[,running_var]>=grid[3] & grid[4]>=d_match2[,running_var],],
             aes(x=running_var,y=y),fill=col2, color=col2, size=2, shape=21,alpha=0.2)+
  geom_point(data=d_match2[d_match2[,running_var]>=grid[2] & grid[3]>=d_match2[,running_var],],
             aes(x=running_var,y=y),fill=col2, color='light grey', size=2, shape=21,alpha=1)+
  geom_vline(aes(xintercept=grid[1]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[2]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[3]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[4]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[5]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[6]),linetype='dashed',col="dark grey",lwd=1) +
  geom_vline(aes(xintercept=grid[8]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[9]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[10]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=grid[11]),linetype='dashed',col="dark grey",lwd=1) + 
  geom_vline(aes(xintercept=0),linetype='dashed',lwd=1.5) +
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="dark grey"), size=1.5)+
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="black"), size=1.5)+
  theme_ipsum_fsc(grid="Y", plot_title_size = 20,
                  axis_text_size = 10,
                  axis_title_size=12,
                  axis_col = col_text)+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  scale_color_identity(name='Potential Outcomes',
                       breaks=c('dark grey','black'),
                       labels = c("Grid","Cutoff"),
                       guide = "legend")+
  #  annotate("pointrage",x=-600,y=41,xend=-600,yend=41,shape=21,size=2,col="dark grey") +
  #  annotate("text",x=13.05,y=8.4,label = "H*", size=5.5, colour=col_text)+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=20, colour=col_text,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12, colour=col_text),
        axis.title.y = element_text(size=20, colour=col_text,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12, colour=col_text),legend.position=c(0.1, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size=15, colour=col_text),
        legend.background = element_rect(fill=alpha("white",0),colour =alpha("white",0)),
        title = element_text(size=25, colour=col_text))+
  ylim(-60,60) + xlim(-600,600)+
  annotate("pointrange",x=-588,y=34,ymin=34,ymax=34,shape=21,size=1,col=col2,fill=col2,alpha=0.2)+
  annotate("text",x=-550,y=34,label = "Matched Obs", size=5, colour=col_text,hjust = 0,family=font_fsc)




ggplot(d_match2,aes(y=y,x=running_var, label=NA), fill=col2)+
#  geom_rect(aes(xmin = grid[2], xmax = grid[9], ymin = -Inf, ymax = Inf),
#            fill = col2, alpha = 0.03) +
  geom_ribbon(aes(x=running_var,ymin = -60, ymax = 60),fill = col2, alpha = 0.3) +
  geom_point(fill=col2, color='light grey', size=2, shape=21)+
  geom_point(data=data2[data2[,tvar]==0,],aes(y=y,x=running_var, label=NA), fill="white",col="slategrey",shape=21,alpha=0.5)+
  geom_vline(aes(xintercept=grid[1]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[2]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[3]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[4]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[5]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[6]),linetype='dashed',col="dark grey",lwd=0.6) +
  geom_vline(aes(xintercept=grid[8]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[9]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[10]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[11]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=0),linetype='dashed',lwd=1.1) +
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="dark grey"), size=1)+
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="black"), size=1)+
  theme_ipsum_fsc(grid="Y", plot_title_size = 20,
                  axis_text_size = 10,
                  axis_title_size=12,
                  axis_col = col_text)+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  scale_color_identity(name='Potential Outcomes',
                       breaks=c('dark grey','black'),
                       labels = c("Grid","Cutoff"),
                       guide = "legend")+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=20, colour=col_text,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12, colour=col_text),
        axis.title.y = element_text(size=20, colour=col_text,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12, colour=col_text),legend.position=c(0.1, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size=15, colour=col_text),
        legend.background = element_rect(fill=alpha("white",0),colour =alpha("white",0)),
        title = element_text(size=25, colour=col_text))+
  ylim(-60,60) + xlim(-600,600)+
  annotate("pointrange",x=-588,y=35,ymin=35,ymax=35,shape=21,size=1,col=col2,fill=col2,alpha=1)+
  annotate("text",x=-550,y=35,label = "Matched Obs", size=5, colour=col_text,hjust = 0,family=font_fsc)+
  annotate("rect",xmin=-600,xmax=-580,ymin=23,ymax=27,col="dark gray",fill=col2,alpha=0.3)+
  annotate("text",x=-550,y=25,label = "Overlap region", size=5, colour=col_text,hjust = 0,family=font_fsc)
  


### Plot 4

plot(data2[data2[,tvar]==0,running_var],data2[data2[,tvar]==0,out_var], pch="o", cex=0.5,col="dark grey",
     ylim=c(-30,30),ann=FALSE)
mtext(side = 1, text = "Running variable", line = 2)
mtext(side = 2, text = "Outcome", line = 2)

points(d_match2[,out_var] ~ d_match2[,running_var],col="deepskyblue3",pch="o",cex=0.6)

polygon(c(rev(eval), eval), c(rev(cil),ciu), 
        col = 'grey80', border = NA)
lines(eval,est,col="darkorchid",lwd=2)
lines(eval, cil, col="black", lty=4)
lines(eval, ciu, col="black", lty=4)

abline(h=est[which(eval==0)],lty=3,col="darkgrey",lwd=2)

for(i in 1:length(grid)){
  abline(v=grid[i],lty=2,col="grey",lwd=1)
}

abline(v=0,col="grey",lty=2,lwd=2)

arrows(-160.1434,23,-240.14340,23,xpd = TRUE,length=0.1)
arrows(96.90316-80,23,96.90316,23,xpd = TRUE,length=0.1)
text(-71,23,bquote(F[j]),xpd=TRUE)

abline(v=optimal_bandwidth[1],col="orange2",lty=6,lwd=2)
abline(v=optimal_bandwidth[2],col="orange2",lty=6,lwd=2)

legend("topright",legend = c("LP regression","95% CI","Point Estimate Y(R=0)","Generalization BW","Matched Obs"),ncol=1,
       lty=c(1,4,3,6,NA),pch=c(NA,NA,NA,NA,"o"),lwd=c(2,1,2,2,NA),
       col=c("darkorchid3","black","darkgrey","orange2","deepskyblue3"),cex=0.8)

data_eval = as.data.frame(cbind(eval,est,cil,ciu))
names(data_eval) = c("running_var","est","cil","ciu")

data3 = data2[data2[,tvar]==0,]

library(gtools)
data3 = smartbind(data3,data_eval)

ggplot(data3,aes(y=y,x=running_var, label=NA), fill="white")+
  #geom_rect(aes(xmin = grid[5], xmax = grid[6], ymin = -Inf, ymax = Inf),
  #          fill = "light grey", alpha = 0.03) +
  geom_ribbon(data=data3[data3$running_var>=grid[5] & data3$running_var<=grid[6],],
              aes(ymin=-30,ymax=30),fill="light grey",alpha=0.5) +
  #geom_point(fill="white", color="slategrey", size=2, shape=21)+
  geom_point(data=d_match2,
             aes(x=running_var,y=y),fill=col2, color=col2, size=2, shape=21,alpha=0.3)+
  geom_vline(aes(xintercept=grid[1]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[2]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[3]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[4]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[5]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[6]),linetype='dashed',col="dark grey",lwd=0.6) +
  geom_vline(aes(xintercept=grid[8]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[9]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[10]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[11]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=0),linetype='dashed',lwd=1) +
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="dark grey"), size=1)+
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="black"), size=1)+
  geom_line(data=data3[!is.na(data3$est),],aes(x=running_var,y=est),col=col_text,lwd=1) +
  geom_ribbon(data=data3[!is.na(data3$est),],aes(x=running_var,ymin=cil,ymax=ciu),fill=col1,alpha=0.3) +
  geom_line(data=data3[!is.na(data3$est),],aes(x=running_var,y=cil), col=col_text, lty=4) +
  geom_line(data=data3[!is.na(data3$est),],aes(x=running_var,y=ciu), col=col_text, lty=4) +
  theme_ipsum_fsc(grid="Y", plot_title_size = 20,
                  axis_text_size = 10,
                  axis_title_size=12,
                  axis_col = col_text)+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  scale_color_identity(name='Potential Outcomes',
                       breaks=c('dark grey','black'),
                       labels = c("Grid","Cutoff"),
                       guide = "legend")+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=20, colour=col_text,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12, colour=col_text),
        axis.title.y = element_text(size=20, colour=col_text,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12, colour=col_text),legend.position=c(0.1, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size=15, colour=col_text),
        legend.background = element_rect(fill=alpha("white",0),colour =alpha("white",0)),
        title = element_text(size=25, colour=col_text))+
  ylim(-30,30) + xlim(-600,600)+
  annotate("pointrange",x=-588,y=18,ymin=18,ymax=18,shape=21,size=1,col=col2,fill=col2,alpha=0.2)+
  annotate("text",x=-550,y=18,label = "Matched Obs", size=5, colour=col_text,hjust = 0,family=font_fsc)+
  annotate("rect",xmin=-600,xmax=-580,ymin=12.3,ymax=13.9,col="dark gray",fill=col1,alpha=0.3)+
  annotate("text",x=-550,y=13.4,label = "Loc. Poly. 95% CI", size=5, colour=col_text,hjust = 0,family=font_fsc) +
  annotate("rect",xmin=-600,xmax=-580,ymin=8,ymax=9.5,col="dark gray",fill="light grey",alpha=0.5)+
  annotate("text",x=-550,y=9,label = "Template 1 region", size=5, colour=col_text,hjust = 0,family=font_fsc)


ggplot(data3,aes(y=y,x=running_var, label=NA), fill="white")+
  #geom_rect(aes(xmin = grid[5], xmax = grid[6], ymin = -Inf, ymax = Inf),
  #          fill = "light grey", alpha = 0.03) +
  geom_ribbon(data=data3[data3$running_var>=grid[5] & data3$running_var<=grid[6],],
              aes(ymin=-30,ymax=30),fill="light grey",alpha=0.5) +
  #geom_point(fill="white", color="slategrey", size=2, shape=21)+
  geom_point(data=d_match2,
             aes(x=running_var,y=y),fill=col2, color=col2, size=2, shape=21,alpha=0.3)+
  geom_vline(aes(xintercept=grid[1]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[2]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[3]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[4]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[5]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[6]),linetype='dashed',col="dark grey",lwd=0.6) +
  geom_vline(aes(xintercept=grid[8]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[9]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[10]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[11]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=optimal_bandwidth[1]),linetype='dashed',col=col1,lwd=1.1) + 
  geom_vline(aes(xintercept=optimal_bandwidth[2]),linetype='dashed',col=col1,lwd=1.1) + 
  geom_vline(aes(xintercept=0),linetype='dashed',lwd=1) +
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="dark grey"), size=1)+
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="black"), size=1)+
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col=col1), size=1)+
  geom_line(data=data3[!is.na(data3$est),],aes(x=running_var,y=est),col=col_text,lwd=1) +
  geom_ribbon(data=data3[!is.na(data3$est),],aes(x=running_var,ymin=cil,ymax=ciu),fill=col1,alpha=0.3) +
  geom_line(data=data3[!is.na(data3$est),],aes(x=running_var,y=cil), col=col_text, lty=4) +
  geom_line(data=data3[!is.na(data3$est),],aes(x=running_var,y=ciu), col=col_text, lty=4) +
  theme_ipsum_fsc(grid="Y", plot_title_size = 20,
                  axis_text_size = 10,
                  axis_title_size=12,
                  axis_col = col_text)+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  scale_color_identity(name='Potential Outcomes',
                       breaks=c('dark grey','black',col1),
                       labels = c("Grid","Cutoff","H1"),
                       guide = "legend")+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=20, colour=col_text,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12, colour=col_text),
        axis.title.y = element_text(size=20, colour=col_text,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12, colour=col_text),legend.position=c(0.1, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size=15, colour=col_text),
        legend.background = element_rect(fill=alpha("white",0),colour =alpha("white",0)),
        title = element_text(size=25, colour=col_text))+
  ylim(-30,30) + xlim(-600,600)+
  annotate("pointrange",x=-588,y=15,ymin=15,ymax=15,shape=21,size=1,col=col2,fill=col2,alpha=0.2)+
  annotate("text",x=-550,y=15,label = "Matched Obs", size=5, colour=col_text,hjust = 0,family=font_fsc)+
  annotate("rect",xmin=-600,xmax=-580,ymin=10,ymax=11.5,col="dark gray",fill=col1,alpha=0.3)+
  annotate("text",x=-550,y=10.5,label = "Loc. Poly. 95% CI", size=5, colour=col_text,hjust = 0,family=font_fsc) +
  annotate("rect",xmin=-600,xmax=-580,ymin=5.2,ymax=6.8,col="dark gray",fill="light grey",alpha=0.5)+
  annotate("text",x=-550,y=6,label = "Template region", size=5, colour=col_text,hjust = 0,family=font_fsc)



## RUN ALL ITERATIONS

random_t = sample(1:nrow(d_match),600,replace=FALSE)
d_match2 = d_match[random_t,]

data_eval_1 = data_eval

data_eval = as.data.frame(cbind(eval,est,cil,ciu))
names(data_eval) = c("running_var","est","cil","ciu")

data_eval_1_trim = data_eval_1[data_eval_1$running_var<max(data_eval$running_var),]
  
data3 = data2[data2[,tvar]==0,]

library(gtools)
#data3 = smartbind(data3,data_eval)
data3 = smartbind(data3,data_eval_1_trim)



ggplot(data3,aes(y=y,x=running_var, label=NA), fill="white")+
  #geom_rect(aes(xmin = bandwidth_opt[1]+15, xmax = bandwidth_opt[2], ymin = -Inf, ymax = Inf),
  #          fill = col3, alpha = 0.003) +
  geom_ribbon(data=data3[data3$running_var>=bandwidth_opt[1]+15 & data3$running_var<=grid[6],],
              aes(ymin=-30,ymax=30),fill="light grey",alpha=0.5) +
 # geom_point(fill="white", color="slategrey", size=2, shape=21)+
  geom_point(data=d_match2,
             aes(x=running_var,y=y),fill=col2, color=col2, size=2, shape=21,alpha=0.2)+
  geom_vline(aes(xintercept=grid[1]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[2]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[3]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[4]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[5]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[6]),linetype='dashed',col="dark grey",lwd=0.6) +
  geom_vline(aes(xintercept=grid[8]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[9]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[10]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=grid[11]),linetype='dashed',col="dark grey",lwd=0.6) + 
  geom_vline(aes(xintercept=optimal_bandwidth[1]+15),linetype='dashed',col=col1,lwd=1.1) + 
  geom_vline(aes(xintercept=optimal_bandwidth[2]),linetype='dashed',col=col1,lwd=1.1) + 
  geom_vline(aes(xintercept=0),linetype='dashed',lwd=1) +
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="dark grey"), size=1)+
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col="black"), size=1)+
  geom_segment(aes(x = -700, y = 0, xend = -701, yend = 0,col=col1), size=1.5)+
  geom_line(data=data3[!is.na(data3$est),],aes(x=running_var,y=est),col=col_text,lwd=1) +
  geom_ribbon(data=data3[!is.na(data3$est),],aes(x=running_var,ymin=cil,ymax=ciu),fill=col1,alpha=0.3) +
  geom_line(data=data3[!is.na(data3$est),],aes(x=running_var,y=cil), col=col_text, lty=4) +
  geom_line(data=data3[!is.na(data3$est),],aes(x=running_var,y=ciu), col=col_text, lty=4) +
  theme_ipsum_fsc(grid="Y", plot_title_size = 20,
                  axis_text_size = 10,
                  axis_title_size=12,
                  axis_col = col_text)+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Pre-intervention")+
  scale_color_identity(name='Potential Outcomes',
                       breaks=c('dark grey','black',col1),
                       labels = c("Grid","Cutoff","H*"),
                       guide = "legend")+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=20, colour=col_text,margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12, colour=col_text),
        axis.title.y = element_text(size=20, colour=col_text,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12, colour=col_text),legend.position=c(0.1, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size=15, colour=col_text),
        legend.background = element_rect(fill=alpha("white",0),colour =alpha("white",0)),
        title = element_text(size=25, colour=col_text))+
  ylim(-30,30) + xlim(-600,600)+
  annotate("pointrange",x=-588,y=15,ymin=15,ymax=15,shape=21,size=1,col=col2,fill=col2,alpha=0.2)+
  annotate("text",x=-550,y=15,label = "Matched Obs", size=5, colour=col_text,hjust = 0,family=font_fsc)+
  annotate("rect",xmin=-600,xmax=-580,ymin=10,ymax=11.5,col="dark gray",fill=col1,alpha=0.3)+
  annotate("text",x=-550,y=10.5,label = "Loc. Poly. 95% CI", size=5, colour=col_text,hjust = 0,family=font_fsc) +
  annotate("rect",xmin=-600,xmax=-580,ymin=5.2,ymax=6.8,col="dark gray",fill="light grey",alpha=0.5)+
  annotate("text",x=-550,y=6,label = "Template region", size=5, colour=col_text,hjust = 0,family=font_fsc)

