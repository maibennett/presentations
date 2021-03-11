# Matching

grd <- function(data, tvar, grid, out_var, cov1, cov2, tols, running_var, threshold, treat, bandwidth,
                size = 1000, nsamples = 10, seed = 100, dir_data, t_max_alloc, match_prop = 0.99,
                significance, drop_tail = 0.01, print=FALSE, name_plot = "plot"){

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
    
    return(list(d_rematch_temp = d_rematch_bef, d_rematch = d_match,
                optimal_bandwidth = optimal_bandwidth, template_n = template_n,
                bias = bias, obw_grid = obw_grid,STOP=STOP))
    }
  }
  
  if(STOP==TRUE){
    print(paste("############################################"))
    print(paste("############################################"))
    print(paste("############################ Error in LPROBUST!"))
    print(paste("############################################"))
    print(paste("############################################"))
    
    return(list(STOP=TRUE))
    
  }
  
}
