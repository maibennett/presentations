#Clear memory
rm(list = ls())

#Clear the console
cat("\014")

# Paper: How far is too far?
# Authors: M. Bennett
# Date: 01-17-19

## Load libraries
library(designmatch)
library(Rglpk)
library(gurobi)
library(xtable)
library(dplyr)
library(readstata13)
library(rdrobust)
library(exact2x2)
library("causaldrf")
library("BayesTree")
library(nprobust)

packages <- c("devtools"
              ,"randomForest"
              ,"rpart" # decision tree
              ,"rpart.plot" # enhanced tree plots
              ,"ROCR"
              ,"Hmisc"
              ,"corrplot"
              ,"texreg"
              ,"glmnet"
              ,"reshape2"
              ,"knitr"
              ,"xtable"
              ,"lars"
              ,"ggplot2"
              ,"matrixStats"
              ,"plyr"
              ,"stargazer")


# Install packages if not installed
list.of.packages <- packages
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

# Load all packages
sapply(packages, require, character.only = TRUE)

#################################################################################################
#################################################################################################

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

dir_data = dirname(rstudioapi::getSourceEditorContext()$path)

#source(paste0(dir_data,"/create_template_func_updated.R"))
source(paste0(dir_data,"/functions/create_template_func_sim3.R"))

load(paste0(dir_data,"/data_clean.Rdata"))

Z = d[d$year==2016 & d$decil_se_f<=9 & d$percapita_f>0,
      c("female_3","edf_8","edm_8","leng_score_11","math_score_11","gpa_score_5",
        "simce_stu","simce_sch","ses_sch_5","rank_score", "gpa_score", "leng_score",
        "math_score","capital", "public_sch","pub_health_2")]


d = d[d$decil_se_f<=9 & d$percapita_f>0,]

# Create variables that are missing:
d_gpa_score_5 <- model.matrix(~ as.factor(gpa_score_5),d)
d_gpa_score_5 <- d_gpa_score_5[,-1]

d_PSUL_11 <- model.matrix(~ as.factor(leng_score_11),d)
d_PSUL_11<-d_PSUL_11[,-1]

d_PSUM_11 <- model.matrix(~ as.factor(math_score_11),d)
d_PSUM_11<-d_PSUM_11[,-1]

names_before = names(d)
d = cbind(d,d_gpa_score_5,d_PSUL_11,d_PSUM_11)
names(d) = c(names_before, "gpa2","gpa3","gpa4","gpa5","PSUL_2","PSUL_3","PSUL_4",
             "PSUL_5","PSUL_6","PSUL_7","PSUL_8","PSUL_9","PSUL_10",
             "PSUM_2","PSUM_3","PSUM_4","PSUM_5","PSUM_6","PSUM_7","PSUM_8",
             "PSUM_9","PSUM_10")

d$treat = as.numeric(d$decil_se_f<=5)
d$after = 0
d$after[d$year==2016]=1

# Drop 2014 for this:

d_pre = d[d$year!=2016,]
d_pre$after = 0
d_pre$after[d_pre$year==2015]=1

d=d[d$year!=2014,]

######################## Estimate covariates that are important:
working = d[d$year==2015,]

covariate.names <- c("female_3","edm_8_lessHS","edm_8_HS","edm_8_lessTec","edm_8_Tec",
                     "edm_8_lessUni","edm_8_Uni","edm_8_miss","edf_8_lessHS","edf_8_HS",
                     "edf_8_lessTec","edf_8_Tec",
                     "edf_8_lessUni","edf_8_Uni","edf_8_miss","PSUL_2","PSUL_3","PSUL_4",
                     "PSUL_5","PSUL_6","PSUL_7","PSUL_8","PSUL_9","PSUL_10",
                     "PSUM_2","PSUM_3","PSUM_4","PSUM_5","PSUM_6","PSUM_7","PSUM_8",
                     "PSUM_9","PSUM_10","gpa2","gpa3","gpa4","gpa5","simce_stu","simce_sch",
                     "ses_sch_5","rank_score", "gpa_score", "leng_score","math_score","capital", 
                     "public_sch","pub_health_2")

# The dependent (outcome) variable is whether the person voted,
# so let's rename "outcome_voted" to Y
names(working)[names(working)=="enrolled_uni"] <- "Y"

# Extract the dependent variable
Y <- working[["Y"]]

# The treatment is whether they received the "your neighbors are voting" letter
names(working)[names(working)=="treat"] <- "W"

# Extract treatment variable & covariates
W <- working[["W"]]
covariates <- working[covariate.names]

# some algorithms require our covariates be scaled
# scale, with default settings, will calculate the mean and standard deviation of the entire vector,
# then "scale" each element by those values by subtracting the mean and dividing by the sd
covariates.scaled <- scale(covariates)
processed.unscaled <- data.frame(Y, W, covariates)
processed.scaled <- data.frame(Y, W, covariates.scaled)

set.seed(44)
smplmain <- sample(nrow(processed.scaled), round(9*nrow(processed.scaled)/10), replace=FALSE)

processed.scaled.train <- processed.scaled[smplmain,]
processed.scaled.test <- processed.scaled[-smplmain,]

y.train <- as.matrix(processed.scaled.train$Y, ncol=1)
y.test <- as.matrix(processed.scaled.test$Y, ncol=1)

# create 45-45-10 sample
smplcausal <- sample(nrow(processed.scaled.train),
                     round(5*nrow(processed.scaled.train)/10), replace=FALSE)
processed.scaled.train.1 <- processed.scaled.train[smplcausal,]
processed.scaled.train.2 <- processed.scaled.train[-smplcausal,]

print(covariate.names)
sumx <- paste(covariate.names, collapse = " + ")  # "X1 + X2 + X3 + ..." for substitution later
interx <- paste(" (",sumx, ")^2", sep="")  # "(X1 + X2 + X3 + ...)^2" for substitution later

# Y ~ X1 + X2 + X3 + ...
linearnotreat <- paste("Y",sumx, sep=" ~ ")
linearnotreat <- as.formula(linearnotreat)
linearnotreat

# Y ~ W + X1 + X2 + X3 + ...
linear <- paste("Y",paste("W",sumx, sep=" + "), sep=" ~ ")
linear <- as.formula(linear)
linear

# Y ~ W * (X1 + X2 + X3 + ...)
# ---> X*Z means include these variables plus the interactions between them
linearhet <- paste("Y", paste("W * (", sumx, ") ", sep=""), sep=" ~ ")
linearhet <- as.formula(linearhet)
linearhet

# LASSO takes in a model.matrix
# First parameter is the model (here we use linear, which we created before)
# Second parameter is the dataframe we want to creaate the matrix from
linear.train <- model.matrix(linearnotreat, processed.scaled.train)[,-1]
linear.test <- model.matrix(linearnotreat, processed.scaled.test)[,-1]
linear.train.1 <- model.matrix(linearnotreat, processed.scaled.train.1)[,-1]
linear.train.2 <- model.matrix(linearnotreat, processed.scaled.train.2)[,-1]


# This is almost the same as the previous part, except that we use alpha in (0, 1)
elastNet.logit <- glmnet(linear.train.1, y.train[smplcausal,],
                         alpha = 0.2, family = 'binomial')

# plot the coefficient paths against the log-lambda values
plot(elastNet.logit, xvar = "lambda", label = FALSE)
grid()

# We can also use the same folds so we can select the optimal value for alpha.
# Above, we set alpha = 0.2
# Use foldid: a vector of values between 1 and nfold identifying what fold
# each observation is in.
foldid <- sample(1:10, size = length(y.train[smplcausal,]), replace = TRUE)
cv0.2.elastNet.logit <- cv.glmnet(linear.train.1, y.train[smplcausal,],
                                  foldid = foldid, alpha = 0.2, parallel=TRUE)
cv0.5.elastNet.logit <- cv.glmnet(linear.train.1, y.train[smplcausal,],
                                  foldid = foldid, alpha = 0.5, parallel=TRUE)
cv0.8.elastNet.logit <- cv.glmnet(linear.train.1, y.train[smplcausal,],
                                  foldid = foldid, alpha = 0.8, parallel=TRUE)

# plot all three MSE's in the same plot to compare
# par(mfrow = c(2,2))
# plot(cv0.2.elastNet.logit); plot(cv0.5.elastNet.logit); plot(cv0.8.elastNet.logit)
plot(log(cv0.8.elastNet.logit$lambda), cv0.8.elastNet.logit$cvm, pch = 19,
     col = "red", xlab = "log(Lambda)", ylab = cv0.2.elastNet.logit$name)
points(log(cv0.5.elastNet.logit$lambda), cv0.5.elastNet.logit$cvm, pch=19, col="grey")
points(log(cv0.2.elastNet.logit$lambda), cv0.2.elastNet.logit$cvm, pch=19, col="blue")
legend("topleft",legend = c("alpha = 0.8", "alpha = 0.5", "alpha = 0.2"),
       pch = 19, col = c("red","grey","blue"))

# We can plot with more values of alpha to choose the best alpha. According to the plot,
# it seems like alpha = 0.8 does the best among the three reported values.

# Using this alpha, we can find the optimal lambda and coefs as in the previous part.
opt0.8.lambda <- cv0.8.elastNet.logit$lambda.1se
coef(cv0.8.elastNet.logit, s = opt0.8.lambda)



# Matching

data = d
tvar = "year2016"
grid = c(quantile(d$inc_pc_adj[d$year2016==1 & d$treat==1],probs=seq(0,1,0.1))[-11],0,
         quantile(d$inc_pc_adj[d$year2016==1 & d$treat==0],probs=seq(0,1,0.1))[-1])
bandwidth = c(max(grid[grid<0]),0)
#cov1 = c("female_3","edm_8_lessHS","edm_8_HS","edm_8_lessTec","edm_8_Tec",
#         "edm_8_lessUni","edm_8_Uni","edm_8_miss","edf_8_lessHS","edf_8_HS",
#         "edf_8_lessTec","edf_8_Tec",
#         "edf_8_lessUni","edf_8_Uni","edf_8_miss","PSUL_2","PSUL_3","PSUL_4",
#         "PSUL_5","PSUL_6","PSUL_7","PSUL_8","PSUL_9","PSUL_10",
#         "PSUM_2","PSUM_3","PSUM_4","PSUM_5","PSUM_6","PSUM_7","PSUM_8",
#         "PSUM_9","PSUM_10","gpa2","gpa3","gpa4","gpa5")
cov1 = c("female_3","edm_8_lessHS","edm_8_HS","edm_8_lessTec","edm_8_Tec",
         "edm_8_lessUni","edm_8_Uni","edm_8_miss","edf_8_lessHS","edf_8_HS",
         "edf_8_lessTec","edf_8_Tec",
         "edf_8_lessUni","edf_8_Uni","edf_8_miss","PSUL_2","PSUL_3","PSUL_4",
         "PSUL_5","PSUL_6","PSUL_7","PSUL_8","PSUL_9","PSUL_10",
         "PSUM_2","PSUM_3","PSUM_4","PSUM_5","PSUM_6","PSUM_7","PSUM_8",
         "PSUM_9","PSUM_10","gpa2","gpa3","gpa4","gpa5")
#cov2 = c("simce_stu","simce_sch","ses_sch_5","rank_score", "gpa_score", "leng_score",
#         "math_score","capital", "public_sch","pub_health_2")
########### variables selected by ML
cov2 = c("simce_stu","ses_sch_5","rank_score", "gpa_score", "leng_score",
         "math_score","capital")
out_var = "enrolled_uni"
#out_var = "applied"
tols = c(0.0001,0.025)
running_var = "inc_pc_adj"
threshold = 0
treat = "treat"
size = 1000
nsamples = 50
seed = 100
dir_data = dir_data
t_max_alloc = 60
significance = 0.1
match_prop = 1
drop_tail=0
#max_bias = 0.01
max_bias = 0.01
print=TRUE
name_plot = "opt_bwv3_enrolled_test"

#### Save original data:
d_input = d

source(paste0(dir_data,"/functions/grd_function_sim3_application_test.R"))

grd_enroll = grd(data = data, tvar = tvar, grid = grid, bandwidth = bandwidth, out_var = out_var, max_bias=max_bias,
                 cov1 = cov1, cov2 = cov2, tols = tols, match_prop = match_prop,
                 running_var = running_var, threshold = threshold, treat = treat, size = size, nsamples = nsamples, 
                 seed = seed, dir_data = dir_data, t_max_alloc = t_max_alloc, significance = significance, drop_tail = drop_tail,
                 print=print, name_plot=name_plot)


#Minimize:
min=FALSE

grid_grd = grd_enroll$optimal_bandwidth

grid0 = list(grid_grd)

if(min==TRUE){
source(paste0(dir_data,"/functions/grd_effects_function3_applicationv2.R"))

grd0 = grd_effect(data = data, tvar = tvar, period=0, grid = grid0, subgroup = FALSE, condition_sg = NULL, 
                  exact=c("female_3","pub_health_2"),
                  out_var=out_var,cov1 = cov1, cov2 = cov2, tols = tols, 
                  running_var = running_var, threshold = threshold, treat = treat, template_grd = grd_enroll$template_n,
                  size = size, nsamples = nsamples, seed = 100, dir_data = dir_data, 
                  t_max_alloc = t_max_alloc) 

}

if(min==FALSE){
  source(paste0(dir_data,"/functions/grd_effects_function3_application.R"))
  
  grd0 = grd_effect(data = data, tvar = tvar, period=0, grid = grid0, subgroup = FALSE, condition_sg = NULL,
                    out_var=out_var,cov1 = cov1, cov2 = cov2, tols = tols, 
                    running_var = running_var, threshold = threshold, treat = treat, template_grd = grd_enroll$template_n,
                    size = size, nsamples = nsamples, seed = 100, dir_data = dir_data, 
                    t_max_alloc = t_max_alloc) 
  
}


if(min==TRUE){
grid_grd = grd_enroll$optimal_bandwidth
grid1 = list(grid_grd)

source(paste0(dir_data,"/functions/grd_effects_function3_applicationv2.R"))

grd1_enrolled = grd_effect(data = data, tvar = tvar, period=1, grid = grid1,exact=c("female_3"),
                           out_var=out_var,cov1 = cov1, cov2 = cov2, tols = tols, 
                           running_var = running_var, threshold = threshold, treat = treat, 
                           template_grd = grd_enroll$template_n,
                           size = size, nsamples = nsamples, seed = 100, dir_data = dir_data, 
                           t_max_alloc = t_max_alloc)
}

if(min==FALSE){

  grid_grd = c(grd_enroll$optimal_bandwidth[1],250)
  
  grid1 = list(grid_grd)
  
  grd1_enrolled = grd_effect(data = data, tvar = tvar, period=1, grid = grid1,
                             out_var=out_var,cov1 = cov1, cov2 = cov2, tols = tols, 
                             running_var = running_var, threshold = threshold, treat = treat, 
                             template_grd = grd_enroll$template_n,
                             size = size, nsamples = nsamples, seed = 100, dir_data = dir_data, 
                             t_max_alloc = t_max_alloc)
}


dat_test = grd1_enrolled$d_rematch[[1]][grd1_enrolled$d_rematch[[1]]$template==0,]

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

#control (161.7) --> USE ORIGINAL 00_optimal_badnwidth.R
mod = lprobust(dat_test$enrolled_uni[dat_test$treat==0 & dat_test$inc_pc_adj<=162.4],
               dat_test$inc_pc_adj[dat_test$treat==0 & dat_test$inc_pc_adj<=162.4])
mod = lprobust(dat_test$enrolled_uni[dat_test$treat==0 & dat_test$inc_pc_adj<=161.7],
               dat_test$inc_pc_adj[dat_test$treat==0 & dat_test$inc_pc_adj<=161.7])
#g2 = nprobust.plot(mod,lcol=col2,CIshade = 0.1,legendGroups = c(bquote(Y^{(1)}~(R)~"|X"))) + 
#  geom_hline(yintercept = 0.476,lty=2,col=col_text,lwd=1.3)
g2 = nprobust.plot(mod,lcol=col2,CIshade = 0.1) + 
  geom_hline(yintercept = 0.476,lty=2,col=col_text,lwd=1.3)
g2 + ylim(0,1) +  #theme_ipsum_fsc(grid="Y", plot_title_size = 20,
                   #               axis_text_size = 10,
                    #              axis_title_size=12,
                    #              axis_col = col_text)+
#  scale_colour_identity(name='Potential Outcomes',breaks=c('Series 1'=col1),
#                        labels = c(bquote(Y^{(1)}~(R)~"|X")), guide="legend")+
  xlab("Running variable") + ylab("Outcome")+ggtitle("Post-intervention")+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(axis.title.x = element_text(size=20, colour=col_text),#margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=12, colour=col_text),
        axis.title.y = element_text(size=20, colour=col_text),#margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y = element_text(size=12, colour=col_text),
        legend.title = element_blank(),legend.position = "none",
        legend.text = element_text(size=15, colour=col_text),
        legend.background = element_rect(fill="white",colour ="white"),
        title = element_text(size=25, colour=col_text))


# For paper

g2 + ylim(0,1) +
  #  scale_colour_identity(name='Potential Outcomes',breaks=c('Series 1'=col1),
  #                        labels = c(bquote(Y^{(1)}~(R)~"|X")), guide="legend")+
  xlab("Adjusted Income per capita") + ylab("Enrollment")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))+
  theme(axis.title.x = element_text(size=12, colour=col_text),#,margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=8, colour=col_text),
        axis.title.y = element_text(size=12, colour=col_text),#,margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size=8, colour=col_text),
        legend.title = element_blank(),legend.position = "none",
        legend.text = element_text(size=15, colour=col_text),
        legend.background = element_rect(fill="white",colour ="white"),
        title = element_text(size=25, colour=col_text))+
  annotate("rect",xmin=125,xmax=130,ymin=0.75,ymax=0.8,col=alpha(col2,0.2),fill=alpha(col2,0.2))+
  annotate("text",x=132,y=0.78,label="95% CI",hjust=0)+
  annotate("segment",x=125,xend=130,y=0.85,yend=0.85,col=col2)+
  annotate("text",x=132,y=0.85,label="Loc. Polynomial",hjust=0)
  
  

#treat
mod = lprobust(dat_test$enrolled_uni[dat_test$treat==1],dat_test$inc_pc_adj[dat_test$treat==1])
nprobust.plot(mod)

summary(mod)
  
diff_method =  0.635 - 0.472

#rd_enrolled = rdrobust(data$enrolled_uni[data$year2016==1],data$inc_pc_adj[data$year2016==1],
#              covs = data[data$year2016==1,c(cov1,cov2)])
rd_enrolled = rdrobust(data$enrolled_uni[data$year2016==1],data$inc_pc_adj[data$year2016==1],
                       covs = Z)

rd_enrolled_plot = rdplot(data$enrolled_uni[data$year2016==1],data$inc_pc_adj[data$year2016==1])

est = list()

beta = grd1_enrolled$results[[1]]

beta0 = grd0$results[[1]]
    
bias2 = mean(beta0$treat-beta0$control)
    
names(beta0) <- c("id","treat_0","control_0")

beta_merged = merge(beta,beta0,by="id",all = FALSE)
    
t = t.test(beta$treat,beta$control,paired=TRUE, alternative = c("two.sided"))
t_dd = t.test(beta_merged$treat-beta_merged$control,beta_merged$treat_0-beta_merged$control_0,
                  paired=TRUE, alternative = c("two.sided"))
  
p = t$estimate
ci.u = t$conf.int[2]
ci.l = t$conf.int[1]
  
est[[1]] = c(p,ci.l,ci.u)
  
p = t_dd$estimate
ci.u = t_dd$conf.int[2]
ci.l = t_dd$conf.int[1]
    
est[[2]]=c(p,ci.l,ci.u)

est_enrolled = est


#### Sensitivity analysis
library(rbounds)

beta_merged$b_i = (0-1)*(0-1)
beta_merged$v_i = -1
beta_merged$s_i = sign(beta_merged$control-beta_merged$treat-beta_merged$b_i*(beta_merged$control_0-beta_merged$treat_0))
beta_merged$d_i = beta_merged$v_i*(beta_merged$control-beta_merged$treat-beta_merged$b_i*(beta_merged$control_0-beta_merged$treat_0))

mcnemar_stat = sum((beta_merged$s_i*beta_merged$v_i+1)/2)

#both pairs have to be discordant:
beta_merged$discordant = 0
beta_merged$discordant[(beta_merged$treat-beta_merged$control)!=0 & (beta_merged$treat_0-beta_merged$control_0)!=0]=1

beta_merged$treated_unit = 0
beta_merged$treated_unit[(beta_merged$treat==1 & beta_merged$control_0==1) | (beta_merged$control==1 & beta_merged$treat_0==1)]=1

beta_merged$J = 0
beta_merged$J[beta_merged$discordant==1 & beta_merged$treated_unit==1] = 1

beta_merged$Di = beta_merged$treat - beta_merged$control - beta_merged$treat_0 + beta_merged$control_0

mcnemar_stat = sum(beta_merged$Di[beta_merged$J==1]/4+1/2)
  
tab = table(beta_merged$discordant,beta_merged$treated_unit)

tab2 = table(beta_merged$s_i,beta_merged$v_i)

s = binarysens(tab2[3],tab2[2],Gamma=2.7,GammaInc = 0.01)

tau_c = sqrt(max(s$bounds$Gamma[s$bounds$`Upper bound`<0.05]))

p_t = tau_c/(tau_c+1)
p_c = 1-p_t


### Test for time invariance:
d_match = grd1_enrolled$d_rematch[[1]]
d_match = d_match[d_match$template==0 & d_match$treat==0,]

eval = c(quantile(d_match[d_match[,running_var]<0,running_var],probs = seq(0,1,0.01))[-length(seq(0,1,0.01))],
         0)

eval = c(0,
         quantile(d_match[d_match[,running_var]>=0,running_var],probs = seq(0,1,0.01))[-1])

lp.model <- lprobust(d_match[,out_var],d_match[,running_var],eval = eval)
nprobust.plot(lp.model)

######################### Angrist & Rokannen (2015)

d_aux_treat = d[d$decil_se4<=5,c("enrolled_uni","inc_pc_adj","female_3","edf_8","edm_8","leng_score_11","math_score_11","gpa_score_5",
                           "simce_stu","simce_sch","ses_sch_5","rank_score", "gpa_score", "leng_score",
                           "math_score","capital", "public_sch","pub_health_2")]

lm_treat = lm(enrolled_uni ~ ., data = d_aux_treat[,-2])
res_treat = lm_treat$residuals
lm_res_treat = lm(res_treat ~ inc_pc_adj, data=d_aux_treat)

######### For Applied

out_var = "applied"

grid_grd = grd_enroll$optimal_bandwidth

grid0 = list(grid_grd)

grd0_app = grd_effect(data = data, tvar = tvar, period=0, grid = grid0,
                  out_var=out_var,cov1 = cov1, cov2 = cov2, tols = tols, 
                  running_var = running_var, threshold = threshold, treat = treat, template_grd = grd_enroll$template_n,
                  size = size, nsamples = nsamples, seed = 100, dir_data = dir_data, 
                  t_max_alloc = t_max_alloc) 


grid_grd = c(grd_enroll$optimal_bandwidth[1],250)

grid1 = list(grid_grd)

grd1_applied = grd_effect(data = data, tvar = tvar, period=1, grid = grid1,
                          out_var=out_var,cov1 = cov1, cov2 = cov2, tols = tols, 
                          running_var = running_var, threshold = threshold, treat = treat, 
                          template_grd = grd_enroll$template_n,
                          size = size, nsamples = nsamples, seed = 100, dir_data = dir_data, 
                          t_max_alloc = t_max_alloc)

est = list()

beta = grd1_applied$results[[1]]

beta0 = grd0_app$results[[1]]
    
bias2 = mean(beta0$treat-beta0$control)
    
names(beta0) <- c("id","treat_0","control_0")
    
beta_merged = merge(beta,beta0,by="id",all = FALSE)

t = t.test(beta$treat,beta$control,paired=TRUE, alternative = c("two.sided"))
t_dd = t.test(beta_merged$treat-beta_merged$control,beta_merged$treat_0-beta_merged$control_0,
                  paired=TRUE, alternative = c("two.sided"))

p = t$estimate
ci.u = t$conf.int[2]
ci.l = t$conf.int[1]
  
est[[1]] = c(p,ci.l,ci.u)
  
p = t_dd$estimate
ci.u = t_dd$conf.int[2]
ci.l = t_dd$conf.int[1]
    
est[[2]]=c(p,ci.l,ci.u)

est_applied = est


#### Sensitivity analysis
library(rbounds)

beta_merged$b_i = (0-1)*(0-1)
beta_merged$v_i = -1
beta_merged$s_i = sign(beta_merged$control-beta_merged$treat-beta_merged$b_i*(beta_merged$control_0-beta_merged$treat_0))
beta_merged$d_i = beta_merged$v_i*(beta_merged$control-beta_merged$treat-beta_merged$b_i*(beta_merged$control_0-beta_merged$treat_0))

mcnemar_stat = sum((beta_merged$s_i*beta_merged$v_i+1)/2)

#both pairs have to be discordant:
beta_merged$discordant = 0
beta_merged$discordant[(beta_merged$treat-beta_merged$control)!=0 & (beta_merged$treat_0-beta_merged$control_0)!=0]=1

beta_merged$treated_unit = 0
beta_merged$treated_unit[(beta_merged$treat==1 & beta_merged$control_0==1) | (beta_merged$control==1 & beta_merged$treat_0==1)]=1

tab = table(beta_merged$discordant,beta_merged$treated_unit)

tab2 = table(beta_merged$s_i,beta_merged$v_i)

s = binarysens(tab2[3],tab2[2],Gamma=2.7,GammaInc = 0.01)

tau_c = sqrt(max(s$bounds$Gamma[s$bounds$`Upper bound`<0.05]))

p_t = tau_c/(tau_c+1)
p_c = 1-p_t


### Test for time invariance:
d_match = grd1_applied$d_rematch[[1]]
d_match = d_match[d_match$template==0 & d_match$treat==0,]

eval = c(quantile(d_match[d_match[,running_var]<0,running_var],probs = seq(0,1,0.01))[-length(seq(0,1,0.01))],
         0)

eval = c(0,
         quantile(d_match[d_match[,running_var]>=0,running_var],probs = seq(0,1,0.01))[-1])

lp.model <- lprobust(d_match[,out_var],d_match[,running_var],eval = eval)
nprobust.plot(lp.model)

### Applied

data = d
tvar = "year2016"
grid = c(quantile(d$inc_pc_adj[d$year2016==1 & d$treat==1],probs=seq(0,1,0.1))[-11],0,
         quantile(d$inc_pc_adj[d$year2016==1 & d$treat==0],probs=seq(0,1,0.1))[-1])
bandwidth = c(max(grid[grid<0]),0)
cov1 = c("female_3","edm_8_lessHS","edm_8_HS","edm_8_lessTec","edm_8_Tec",
         "edm_8_lessUni","edm_8_Uni","edm_8_miss","edf_8_lessHS","edf_8_HS",
         "edf_8_lessTec","edf_8_Tec",
         "edf_8_lessUni","edf_8_Uni","edf_8_miss","PSUL_2","PSUL_3","PSUL_4",
         "PSUL_5","PSUL_6","PSUL_7","PSUL_8","PSUL_9","PSUL_10",
         "PSUM_2","PSUM_3","PSUM_4","PSUM_5","PSUM_6","PSUM_7","PSUM_8",
         "PSUM_9","PSUM_10","gpa2","gpa3","gpa4","gpa5")
cov2 = c("simce_stu","simce_sch","ses_sch_5","rank_score", "gpa_score", "leng_score",
         "math_score","capital", "public_sch","pub_health_2")
#out_var = "enrolled_uni"
out_var = "applied"
tols = c(0.0001,0.025)
running_var = "inc_pc_adj"
threshold = 0
treat = "treat"
size = 1000
nsamples = 50
seed = 100
dir_data = dir_data
t_max_alloc = 60
significance = 0.1
match_prop = 1
drop_tail=0
print=TRUE
name_plot = "opt_bwv3_applied"

#### Save original data:
d_input = d

source(paste0(dir_data,"/functions/grd_function_sim3_application.R"))

grd_applied = grd(data = data, tvar = tvar, grid = grid, bandwidth = bandwidth, out_var = out_var, 
                  cov1 = cov1, cov2 = cov2, tols = tols, match_prop = match_prop,
                  running_var = running_var, threshold = threshold, treat = treat, size = size, nsamples = nsamples, 
                  seed = seed, dir_data = dir_data, t_max_alloc = t_max_alloc, significance = significance, drop_tail = drop_tail,
                  print=print, name_plot=name_plot)


source(paste0(dir_data,"/functions/grd_effects_function3_application.R"))

grid_grd = c(-500.26,250)

grid0 = list(grid_grd)

grd0_app2 = grd_effect(data = data, tvar = tvar, period=0, grid = grid0,
                  out_var=out_var,cov1 = cov1, cov2 = cov2, tols = tols, 
                  running_var = running_var, threshold = threshold, treat = treat, template_grd = grd_applied$template_n,
                  size = size, nsamples = nsamples, seed = 100, dir_data = dir_data, 
                  t_max_alloc = t_max_alloc, new_template=TRUE) 


grid_narrow = c(max(grid[grid<0]),min(grid[grid>0]))
grid_grd = c(-500.26,250)

grid1 = list(grid_narrow, grid_grd)

grd1_applied = grd_effect(data = data, tvar = tvar, period=1, grid = grid1,
                          out_var=out_var,cov1 = cov1, cov2 = cov2, tols = tols, 
                          running_var = running_var, threshold = threshold, treat = treat, 
                          template_grd = grd0_app$template,
                          size = size, nsamples = nsamples, seed = 100, dir_data = dir_data, 
                          t_max_alloc = t_max_alloc)


rd_applied = rdrobust(data$applied[data$year2016==1],data$inc_pc_adj[data$year2016==1],
                      covs = Z)

rd_applied_plot = rdplot(data$applied[data$year2016==1],data$inc_pc_adj[data$year2016==1])

head(rd_applied_plot$vars_poly[rd_applied_plot$vars_poly$rdplot_x>=0,])


beta = grd1_applied$results[[2]]
  
beta0 = grd0_app$results[[1]]
    
bias2 = mean(beta0$treat-beta0$control)
    
names(beta0) <- c("id","treat_0","control_0")
    
beta_merged = merge(beta,beta0,by="id",all = FALSE)
    
t = t.test(beta$treat,beta$control,paired=TRUE, alternative = c("two.sided"))
t_dd = t.test(beta_merged$treat-beta_merged$control,beta_merged$treat_0-beta_merged$control_0,
              paired=TRUE, alternative = c("two.sided"))

p = t$estimate
ci.u = t$conf.int[2]
ci.l = t$conf.int[1]
  
est[[1]] = c(p,ci.l,ci.u)
  
p = t_dd$estimate
ci.u = t_dd$conf.int[2]
ci.l = t_dd$conf.int[1]
    
est[[2]]=c(p,ci.l,ci.u)

est_applied2 = est


########################################################################################################
### Pre-intervention DD
########################################################################################################

grid_grd = c(grd_enroll$optimal_bandwidth[1],100)

grid1 = list(grid_grd)

grd0_enrolled_pre = grd_effect(data = d_pre, tvar = "year2015", period=0, grid = grid1,
                               out_var="enrolled_uni",cov1 = cov1, cov2 = cov2, tols = tols, 
                               running_var = running_var, threshold = threshold, treat = treat, 
                               template_grd = grd_enroll$template_n,
                               size = size, nsamples = nsamples, seed = 100, dir_data = dir_data, 
                               t_max_alloc = t_max_alloc)


grd0_applied_pre = grd_effect(data = d_pre, tvar = "year2015", period=0, grid = grid1,
                              out_var="applied",cov1 = cov1, cov2 = cov2, tols = tols, 
                              running_var = running_var, threshold = threshold, treat = treat, 
                              template_grd = grd_enroll$template_n,
                              size = size, nsamples = nsamples, seed = 100, dir_data = dir_data, 
                              t_max_alloc = t_max_alloc)


est = list()

for(s in 1:length(grid1)){
  
  print(s)
  
  beta = grd0$results[[s]]
  
  beta0 = grd0_enrolled_pre$results[[1]]
  
  names(beta0) <- c("id","treat_0","control_0")
  
  beta_merged = merge(beta,beta0,by="id",all = FALSE)
  
  t = t.test(beta$treat,beta$control,paired=TRUE, alternative = c("two.sided"))
  t_dd = t.test(beta_merged$treat-beta_merged$control,beta_merged$treat_0-beta_merged$control_0,
                paired=TRUE, alternative = c("two.sided"))
  
  p = t$estimate
  ci.u = t$conf.int[2]
  ci.l = t$conf.int[1]
  
  est[[s]] = c(p,ci.l,ci.u)
  
  p = t_dd$estimate
  ci.u = t_dd$conf.int[2]
  ci.l = t_dd$conf.int[1]
  
  est[[2]]=c(p,ci.l,ci.u)
}

est_enrolled_pre = est


est = list()

for(s in 1:length(grid1)){
  
  print(s)
  
  beta = grd0_app$results[[s]]
  
  beta0 = grd0_applied_pre$results[[1]]
  
  names(beta0) <- c("id","treat_0","control_0")
  
  beta_merged = merge(beta,beta0,by="id",all = FALSE)
  
  t = t.test(beta$treat,beta$control,paired=TRUE, alternative = c("two.sided"))
  t_dd = t.test(beta_merged$treat-beta_merged$control,beta_merged$treat_0-beta_merged$control_0,
                paired=TRUE, alternative = c("two.sided"))
  
  p = t$estimate
  ci.u = t$conf.int[2]
  ci.l = t$conf.int[1]
  
  est[[s]] = c(p,ci.l,ci.u)
  
  p = t_dd$estimate
  ci.u = t_dd$conf.int[2]
  ci.l = t_dd$conf.int[1]
  
  est[[2]]=c(p,ci.l,ci.u)
}

est_applied_pre = est


### Plot results:
years = c(2014.95,2015.05,2015.95,2016.05)
point_est = c(est_applied_pre[[2]][1],est_enrolled_pre[[2]][1],est_applied[[2]][1],est_enrolled[[2]][1])
cil_est = c(est_applied_pre[[2]][2],est_enrolled_pre[[2]][2],est_applied[[2]][2],est_enrolled[[2]][2])
ciu_est = c(est_applied_pre[[2]][3],est_enrolled_pre[[2]][3],est_applied[[2]][3],est_enrolled[[2]][3])

plot(years,
     point_est,
     xlab = "Years",
     ylab = "Difference in outcome",
     pch = "", ylim=c(-0.2,0.2),xaxt="n",xlim=c(2014.5,2016.5))
axis(1,at=c(2015,2016))
points(years[c(1,3)],point_est[c(1,3)],pch=21,col="deepskyblue3",cex=1.5,lwd=2)
points(years[c(2,4)],point_est[c(2,4)],pch=24,col="darkorchid3",cex=1.5,lwd=2)
abline(h=0, lty = 1, col = "black")

# significance % pointwise standard errors
for(i in 1:4){
segments(years[i], cil_est[i],years[i],ciu_est[i], lty = 1, col = 'dark grey',lwd=1.5)
}

legend('bottomright',
       c("Diff T-C Applications","Diff T-C Enrollment", "95% Confidence Intervals"),
       lty = c(NA,NA,1),
       pch=c(21,24,NA),
       bty = 'Outcome',
       cex = 1,
       col = c("deepskyblue3","darkorchid3", "dark grey"),
       lwd = 2,
       bg = "white")



### Plot results:

d_obw = d[(d$inc_pc_adj>=grd_enroll$optimal_bandwidth[1] & d$inc_pc_adj<=grd_enroll$optimal_bandwidth[2]),]

t_all_bef_enrolled = t.test(d_obw$enrolled_uni[d_obw$year==2015 & d_obw$treat==1],
                            d_obw$enrolled_uni[d_obw$year==2015 & d_obw$treat==0],conf.int = TRUE)

t_all_af_enrolled = t.test(d_obw$enrolled_uni[d_obw$year==2016 & d_obw$treat==1],
                           d_obw$enrolled_uni[d_obw$year==2016 & d_obw$treat==0],conf.int = TRUE)

t_all_bef_applied = t.test(d_obw$applied[d_obw$year==2015 & d_obw$treat==1],
                           d_obw$applied[d_obw$year==2015 & d_obw$treat==0],conf.int = TRUE)

t_all_af_applied = t.test(d_obw$applied[d_obw$year==2016 & d_obw$treat==1],
                          d_obw$applied[d_obw$year==2016 & d_obw$treat==0],conf.int = TRUE)

point_est = c(t_all_bef_applied$estimate[1]-t_all_bef_applied$estimate[2],
              t_all_bef_enrolled$estimate[1]-t_all_bef_enrolled$estimate[2],
              t_all_af_applied$estimate[1]-t_all_af_applied$estimate[2],
              t_all_af_enrolled$estimate[1]-t_all_af_enrolled$estimate[2])

cil_est = c(t_all_bef_applied$conf.int[1],
            t_all_bef_enrolled$conf.int[1],
            t_all_af_applied$conf.int[1],
            t_all_af_enrolled$conf.int[1])
ciu_est = c(t_all_bef_applied$conf.int[2],
            t_all_bef_enrolled$conf.int[2],
            t_all_af_applied$conf.int[2],
            t_all_af_enrolled$conf.int[2])

plot(years,
     point_est,
     xlab = "Years",
     ylab = "Difference in outcome",
     pch = "", ylim=c(-0.2,0.2),xaxt="n",xlim=c(2014.5,2016.5))
axis(1,at=c(2015,2016))
points(years[c(1,3)],point_est[c(1,3)],pch=21,col="deepskyblue3",cex=1.5,lwd=2)
points(years[c(2,4)],point_est[c(2,4)],pch=24,col="darkorchid3",cex=1.5,lwd=2)
abline(h=0, lty = 1, col = "black")

# significance % pointwise standard errors
for(i in 1:4){
  segments(years[i], cil_est[i],years[i],ciu_est[i], lty = 1, col = 'dark grey',lwd=1.5)
}

legend('bottomright',
       c("Diff T-C Applications","Diff T-C Enrollment", "95% Confidence Intervals"),
       lty = c(NA,NA,1),
       pch=c(21,24,NA),
       bty = 'Outcome',
       cex = 1,
       col = c("deepskyblue3","darkorchid3", "dark grey"),
       lwd = 2,
       bg = "white")


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


ggplot(aes(y=years,x=point_est, label=NA), fill="white")+geom_point(fill="white", color="slategrey", size=2, shape=21)+
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
        legend.background = element_rect(fill="white",colour ="white"),
        title = element_text(size=25, colour=col_text))+
  ylim(-60,60) + xlim(-600,600)+
  annotate("pointrange",x=-588,y=40,ymin=40,ymax=40,shape=21,size=1,col="slategrey",fill="white",alpha=1)+
  annotate("text",x=-550,y=40,label = "Observations", size=5, colour=col_text,hjust = 0,family=font_fsc)


###### Check mean balance!

## Comparison between treated observations in the narrow bandwidth and the generalization bandwidth

d$narrow_bandwidth=0
d$narrow_bandwidth[d$inc_pc_adj>=-70 & d$inc_pc_adj<=70 & d$year2016==1] = 1

d$grd_bandwidth=0
d$grd_bandwidth[d$inc_pc_adj>=-500.9 & d$inc_pc_adj<=709 & d$year2016==1] = 1

d_aux = d[d$year2016==1,]

names_mom <- c("inc_pc","female_3","edm","edf","leng_score","math_score","gpa_score","rank_score",
               "simce_stu","simce_sch","ses_sch_5","capital", "public_sch","pub_health_2")

meantab1 = meantab(d_aux[d_aux$narrow_bandwidth==1,names_mom],d_aux$treat[d_aux$narrow_bandwidth==1],
                   which(d_aux$treat[d_aux$narrow_bandwidth==1]==1),
                   which(d_aux$treat[d_aux$narrow_bandwidth==1]==0))[,-c(1,2,3,7)]

row.names(meantab1) = c("Income per capita","Female","Mother's education (years)","Father's education (years)",
                        "Language PSU score", "Math PSU score","GPA score", "Ranking score",
                        "SIMCE 10th grade (student)","SIMCE 10th grade (school)",
                        "SES group school","Lives in Metropolitan region", "Public school",
                        "Public health insurance")


meantab2 = meantab(d_aux[d_aux$grd_bandwidth==1,names_mom],d_aux$treat[d_aux$grd_bandwidth==1],
                   which(d_aux$treat[d_aux$grd_bandwidth==1]==1),
                   which(d_aux$treat[d_aux$grd_bandwidth==1]==0))[,-c(1,2,3,7)]

row.names(meantab2) = c("Income per capita","Female","Mother's education (years)","Father's education (years)",
                        "Language PSU score", "Math PSU score","GPA score", "Ranking score",
                        "SIMCE 10th grade (student)","SIMCE 10th grade (school)",
                        "SES group school","Lives in Metropolitan region", "Public school",
                        "Public health insurance")

star = NA
# differences:
for(i in 1:length(names_mom)){
  x_narrow = d_aux[d_aux$narrow_bandwidth==1,names_mom[i]]
  x_grd = d_aux[d_aux$grd_bandwidth==1,names_mom[i]]
  
  t = t.test(x_narrow, x_grd, alternative = "two.sided", conf.level = 0.95)
  
  if(t$p.value<=0.01){
    star = c(star, "***")
  }
  if(t$p.value>0.01 & t$p.value<=0.05){
    star = c(star, "**")
  }
  
  if(t$p.value>0.05 & t$p.value<=0.1){
    star = c(star, "*")
  }
  if(t$p.value>0.1){
    star = c(star,"")
  }
}

star = star[-1]

tab_diff = cbind(meantab1[,1],meantab2[,1],paste0(round(meantab2[,1]-meantab1[,1],2),star))

#Pre-intervention:

d_enroll = grd0$d_rematch[[1]]

d_enroll = d_enroll[d_enroll$template==0,]

mean_bal_2015 = meantab(d_enroll[,cov2],d_enroll$treat,which(d_enroll$treat==1),which(d_enroll$treat==0))[,4:6]


d_enroll_16 = grd1_enrolled$d_rematch[[1]]

d_enroll_16 = d_enroll_16[d_enroll_16$template==0,]

mean_bal_2016 = meantab(d_enroll_16[,cov2],d_enroll$treat,which(d_enroll_16$treat==1),which(d_enroll_16$treat==0))[,4:6]


xtable(cbind(mean_bal_2015, mean_bal_2016))


finebal = rbind(cbind(table(d_enroll$female_3,d_enroll$treat),
                      table(d_enroll_16$female_3,d_enroll_16$treat)),
                cbind(table(d_enroll$edf_8,d_enroll$treat),
                      table(d_enroll_16$edf_8,d_enroll_16$treat)),
                cbind(table(d_enroll$edm_8,d_enroll$treat),
                      table(d_enroll_16$edm_8,d_enroll_16$treat)),
                cbind(table(d_enroll$leng_score_11,d_enroll$treat),
                      table(d_enroll_16$leng_score_11,d_enroll_16$treat)))


xtable(finebal)

