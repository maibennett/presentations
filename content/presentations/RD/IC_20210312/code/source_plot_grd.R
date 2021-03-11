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
library(Smisc)

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

sim = 33 # This is for the first plots (saved in example_grd.Rdata)

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
dir_data = ""
t_max_alloc = 30
significance = 0.05
match_prop = 0.99
if(sample_size=="small"){match_prop = 0.95}
drop_tail=0
if(sample_size=="small"){drop_tail = 0}
print = FALSE
if(sim==1){print=TRUE}
name_plot = paste0("BW_corr_",correlation,"_size_",sample_size,"_relyr_",rel_y_r,"_unobs_",unobs)


############ grd_function_sim3.R: ONE ITERATION

source("https://raw.githubusercontent.com/maibennett/presentations/main/content/presentations/RD/IC_20210312/code/create_template_func_sim3.R")

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

#### SAvE ENVIRONMENT FOR THIS

#### Now run ALL ITERATIONS

source("https://raw.githubusercontent.com/maibennett/presentations/main/content/presentations/RD/IC_20210312/code/grd_function_sim3.R")

source("https://raw.githubusercontent.com/maibennett/presentations/main/content/presentations/RD/IC_20210312/code/create_template_func_sim3.R")

bandwidth_aux = bandwidth
bandwidth_opt = bandwidth

grid_t = grid[grid<0]
grid_t = grid_t[order(grid[grid<0],decreasing = TRUE)]

optimal_bandwidth = NA

STOP = FALSE
OPT = FALSE
at_least_one = FALSE
J = 1

while(OPT==FALSE & STOP==FALSE){
  
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
  J = J+1
}

if(OPT==TRUE | (OPT==FALSE & at_least_one==TRUE)){
  
  data = as.data.frame(data)
  
  data$id_grd = seq(1,nrow(data),1)
  
  data$template = 0
  
  sig_dif = NA
  
  ct = create_template(data = data, id = "id_grd", running_var = running_var,
                       cov1 = cov1, cov2 = cov2,
                       bandwidth = bandwidth_opt, size = size, 
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
  ecdf_r = quantile(d_match[,"running_var"], probs = c(drop_tail/2,(1-drop_tail/2)))
  low_r = ecdf_r[1]
  up_r = ecdf_r[2]
  
  #d_match = d_match[d_match[,running_var]>=low_r & d_match[,running_var]<=up_r,]
  
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
      bias_u = est[eval>=0][nl[1]]
      obw_grid1_u = gri[max(which(gri/opt_bw1_u<1))]
    }
    
    bias = bias_l - bias_u
    
    optimal_bandwidth = c(max(c(opt_bw1_l,opt_bw2_l)),min(c(opt_bw1_u,opt_bw2_u)))
    obw_grid = c(max(c(obw_grid1_l,obw_grid2_l)),min(c(obw_grid1_u,obw_grid2_u)))
    template_n = template
    
    if(print==TRUE){pdf(paste0(dir_data,"/sim/plots/",name_plot,".pdf"),width=11, height=8.5)}
    
    plot(d_match[,out_var] ~ d_match[,running_var], pch=16, cex=0.6,col="light grey",
         ylim=c(-20,20),xlab="Running Variable",ylab="Outcome", xlim=c(low_r,up_r))
    polygon(c(rev(eval), eval), c(rev(cil),ciu), 
            col = 'grey80', border = NA)
    lines(eval,est,col="blue",lwd=2)
    lines(eval, cil, col="blue", lty=2)
    lines(eval, ciu, col="blue", lty=2)
    
    abline(h=est[which(eval==0)],lty=2,col="black")
    
    for(i in 1:length(grid)){
      abline(v=grid[i],lty=2,col="darkgrey",lwd=1.5)
    }
    
    abline(v=optimal_bandwidth[1],col="red",lty=2)
    abline(v=optimal_bandwidth[2],col="red",lty=2)
    legend("topleft",legend = c("LP regression","95% CI","Point Estimate Y(R=0)","Generalization BW","Grid"),ncol=1,
           lty=c(1,2,2,2,2),lwd=c(2,1,1,1,1.5),col=c("blue","blue","black","red","darkgrey"),cex=0.8)
    
    if(print==TRUE){dev.off()}
    
    if(print==TRUE){pdf(paste0(dir_data,"/sim/plots/",name_plot,"_all.pdf"),width=11, height=8.5)}
    
    plot(data[data[,tvar]==0,running_var],data[data[,tvar]==0,out_var], pch=".", cex=0.8,col="dark grey",
         ylim=c(-60,60),xlab="Running Variable",ylab="Outcome")
    
    points(d_match[,out_var] ~ d_match[,running_var],col="forest green",pch="o",cex=0.4)
    polygon(c(rev(eval), eval), c(rev(cil),ciu), 
            col = 'grey80', border = NA)
    lines(eval,est,col="blue",lwd=2)
    lines(eval, cil, col="blue", lty=2)
    lines(eval, ciu, col="blue", lty=2)
    
    abline(h=est[which(eval==0)],lty=3,col="grey")
    
    for(i in 1:length(grid)){
      abline(v=grid[i],lty=2,col="black",lwd=1.5)
    }
    
    abline(v=0,col="darkorchid3",lty=2,lwd=2)
    
    abline(v=optimal_bandwidth[1],col="red",lty=2)
    abline(v=optimal_bandwidth[2],col="red",lty=2)
    
    legend("topleft",legend = c("LP regression","95% CI","Point Estimate Y(R=0)","Generalization BW","Grid","Matched Obs"),ncol=1,
           lty=c(1,2,2,2,2,NA),pch=c(NA,NA,NA,NA,NA,"o"),lwd=c(2,1,1,1,1.5,NA),col=c("blue","blue","black","red","darkgrey","forest green"),cex=0.8)
    
    if(print==TRUE){dev.off()}
    
    print("#############################################")
    print("#############################################")
    print(paste("############ Sim:",sim))
    print(paste("############ Bandwidth:",round(optimal_bandwidth[1],2),"-",round(optimal_bandwidth[2],2)))
    print("#############################################")
    print("#############################################")
    
    d_rematch_temp = d_rematch_bef 
    d_rematch = d_match
    optimal_bandwidth = optimal_bandwidth
    template_n = template_n
    bias = bias 
    obw_grid = obw_grid
    STOP=STOP
  }
}

if(STOP==TRUE){
  print(paste("############################################"))
  print(paste("############################################"))
  print(paste("############################ Error in LPROBUST!"))
  print(paste("############################################"))
  print(paste("############################################"))
  
  STOP=TRUE
  
}