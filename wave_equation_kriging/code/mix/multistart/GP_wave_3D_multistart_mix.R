### clear workspace
rm(list = ls())

### indicate your workspace and library directory if required (e.g. if code running on server)
### you need to put the correct datasets in the same folder as this file
str_dir <- "/your_work_directory"
str_lib <- "/your_library_directory"

setwd(str_dir)

### uncomment to install the packages in the library directory
#install.packages("remotes",lib=str_lib)
#install.packages("installr",lib=str_lib)
#install.packages("withr",lib=str_lib)
#install.packages("lattice",lib=str_lib)
#install.packages("testthat",lib=str_lib)
#install.packages("nloptr",lib=str_lib)
#install.packages("Rcpp",lib=str_lib)
#install.packages("digest",lib=str_lib)
#install.packages("kergp",lib=str_lib)
#install.packages("foreach",lib=str_lib)
#install.packages("doParallel",lib=str_lib)
#install.packages("lhs",lib=str_lib)
#install.packages("pracma",lib=str_lib)

### import libraries
library("pracma",lib.loc=str_lib)
library("lattice",lib.loc=str_lib)
library("nloptr",lib.loc=str_lib)
library("testthat",lib.loc=str_lib)
library("digest",lib.loc=str_lib)
library("Rcpp",lib.loc=str_lib)
library("kergp",lib.loc=str_lib)
library("foreach",lib.loc=str_lib)
library("doParallel",lib.loc=str_lib)
library("lhs",lib.loc=str_lib)

### import dependencies
source('kernWaveMaternFun_pos_spd.R')  # contains WIGPR covariance function
source('selectPtsToPredict.R')         # small function which discards the informationless points

### Set booleans for the script
use_cpp <- 1

### for reproducibility of latin hypercube
set.seed(1)

### new grid parameters 
n_exp <- 1 # only use 1 experience for multistart procedure
n_multistart <- 100
list_sensors <- c(1,2,3,4,5,6,8,10,15,20,25,30)
len_sensors <- length(list_sensors)

n_space <- 230
n_t <- 2
x_rng <- 1
t_rng <- 1e-7
dt <- t_rng/(n_t-1)

### Build manual kernels for the wave equation
print('Building manual kernels for the wave equation...')

# we use a log scale for the characteristic lenghscales: mu = log(rho)
# order of the hyperparameters : FIRST POSITION then SPEED
parLower_h <- c(ru_x = 0, ru_y = 0, ru_z = 0, Ru = 0.05, mu_u = log(0.01), var_u = 0.1, rv_x = 0, rv_y = 0, rv_z = 0, Rv = 0.05, mu_v = log(0.01), var_v = 0.1, speed = 0.2)
parUpper_h <- c(ru_x = 1, ru_y = 1, ru_z = 1, Ru = 0.4,  mu_u = log(1),    var_u = 5  , rv_x = 1, rv_y = 1, rv_z = 1, Rv = 0.4,  mu_v = log(0.1),    var_v = 5,   speed = 0.8)
par_default <- c(ru_x = 0.65, ru_y = 0.3, ru_z = 0.5, Ru = 0.3, mu_u = log(0.2), var_u = 3, rv_x = 0.3, rv_y = 0.6, rv_z = 0.7, Rv = 0.15, mu_v = log(0.015), var_v = 3, speed = 0.5)

parNames_h <- c("ru_x", "ru_y", "ru_z", "Ru", "log_range_u", "var_u","rv_x", "rv_y", "rv_z", "Rv", "log_range_v", "var_v", "speed")
n_par <- length(par_default)
varNoise_upper <- 2*1e-2
varNoise_lower <- 1e-8

### Define the wave equation kernel for radially symmetric initial position and speed
{
  inputNames_h <- c("x", "y", "z", "t")
  
  covWaveMaternCPP <- covMan(kernel = cppWave3D_pos_spd,
                             acceptMatrix = FALSE,
                             hasGrad = FALSE,
                             label = "non stationary compactly supported wave kernel with cpp implementation",
                             d = 4,
                             inputs = inputNames_h,
                             parLower = parLower_h,
                             parUpper = parUpper_h,
                             parNames = parNames_h,
                             par = par_default)
}

### Build a new regular grid for predictions
{
  myGrid_u <- expand.grid(x = seq(-0, x_rng, length.out = n_space),
                          y = seq(-0, x_rng, length.out = n_space),
                          z = par_default[3], 
                          t = seq(0, t_rng, length.out = n_t), 
                          KEEP.OUT.ATTRS = FALSE)
  myGrid_v <- expand.grid(x = seq(-0, x_rng, length.out = n_space),
                          y = seq(-0, x_rng, length.out = n_space),
                          z = par_default[9], 
                          t = seq(0, t_rng, length.out = n_t), 
                          KEEP.OUT.ATTRS = FALSE)
}

for (id_exp in 1:n_exp){
  
  print(sprintf('exp %s out of %s',num2str(id_exp,fmt=0),num2str(n_exp,fmt=0)))
  
  ### Read data
  wave_full <- read.csv(sprintf("exp_%s_data_set_noisy.csv",num2str(id_exp,fmt=0)), sep = ",", header = FALSE) #"data_set.csv"
  names(wave_full) <- c(inputNames_h, "u")
  nb_t_per_sensor <- 75
  
  ### Define data arrays for storing hyperparam estimation results
  {
    hyp_mat <- matrix(0,len_sensors,14) # format : hyperparameters + estimated noise
    dist_mat <-matrix(0,len_sensors,5) # format : see hyperparameters
    n_err <- matrix(0,len_sensors,1)
    time_hyp <- matrix(0,len_sensors)
    time_err <- matrix(0,len_sensors)
  }
  
  for (id_sensors in 1:len_sensors){
    
    nb_sensors_used <- list_sensors[id_sensors]
    
    print(sprintf('exp %s : nb sensors used : %s',num2str(id_exp,fmt=0),num2str(nb_sensors_used,fmt=0)))
    wave <- wave_full[1:(nb_sensors_used*nb_t_per_sensor),]
    wave_pos <- wave[,1:4]
    n_data <- nrow(wave_pos)
    
    estim_hyp <- TRUE
    var_compact <- 1e-5

    if (estim_hyp){
      
      print(sprintf('exp %s : estimating hyperparameters...',num2str(id_exp,fmt=0)))
      tic()
      
      # generate latin hypercube for multistart procedure
      lh <- maximinLHS(n_multistart, n_par+1, dup=2) # n_par+1 for noise!!
      lh_par <- lh
      for (id_par in 1:n_par){
        lh_par[,id_par] <- parLower_h[id_par] + (parUpper_h[id_par] - parLower_h[id_par]) * lh_par[,id_par] 
      }
      lh_par[,n_par+1] <- varNoise_lower + (varNoise_upper - varNoise_lower)*lh_par[,n_par+1]
      
      npar2 <- n_par+2
      hyp_opt_reshape <- matrix(nrow = n_multistart,ncol = npar2)
      for (id_opt in 1:n_multistart){
        
        parCovInitial <- lh_par[id_opt,1:n_par]
        var_ini <- lh_par[id_opt,n_par+1]
        myGP_dummy <- my_gp(formula = u ~ 1, data = wave, estim = TRUE,
                            cov = covWaveMaternCPP, varNoiseIni = var_ini, varNoiseLower = varNoise_lower, varNoiseUpper = varNoise_upper, parCovIni = parCovInitial, multistart = 1, beta = 0)
        
        hyp_and_loglik <- c(coef(myGP_dummy$covariance),myGP_dummy$varNoise,myGP_dummy$MLE$logLik)
        write.table(hyp_and_loglik, file = sprintf("exp_%s_capt_%s_multi_%s_hyp_and_loglik.csv", 
                                                   num2str(id_exp,fmt=0),num2str(list_sensors[id_sensors],fmt=0), num2str(id_opt,fmt=0)),row.names = FALSE)
        hyp_opt_reshape[id_opt,] <- hyp_and_loglik
      }
      
      colnames(hyp_opt_reshape) <- c(parNames_h,'varNoise','logLik')
      id_best_hyp <- which.max(hyp_opt_reshape[,npar2])
      max_loglik <- max(hyp_opt_reshape[,npar2])
      best_hyp <- hyp_opt_reshape[id_best_hyp,]
      best_noise <- best_hyp[npar2-1]
      coef(covWaveMaternCPP) <- best_hyp[1:n_par]
      
      myGP_compact_cpp <- gp(formula = u ~ 1, data = wave, estim = FALSE, 
                             cov = covWaveMaternCPP, varNoise = best_noise, parCovIni = parCovInitial, beta = 0)
      
      time_hyp[id_sensors] <- as.numeric(toc())
    } 
    
    hyp_h <- coef(myGP_compact_cpp$covariance)
    # for cos
    dist_x0_u <- sqrt(sum((hyp_h[1:3] - par_default[1:3])^2))
    dist_Ru <- abs(hyp_h[4]-par_default[4])
    dist_x0_v <- sqrt(sum((hyp_h[7:9] - par_default[7:9])^2))
    dist_Rv <- abs(hyp_h[10]-par_default[10])
    dist_c <- abs(hyp_h[13]-par_default[13])
    
    dist_h <-c(dist_x0_u,dist_Ru,dist_x0_v,dist_Rv,dist_c)
    
    # Store estimated hyperparameters
    hyp_mat[id_sensors,] <- c(hyp_h,myGP_compact_cpp$varNoise)
    dist_mat[id_sensors,] <- dist_h
    
    #print(hyp_h)
    #print('Estimated noise...')
    #print(myGP_compact_cpp$varNoise)
    #print(dist_h)
    
    # Store estimated hyperparameters
    hyp_mat[id_sensors,] <- c(hyp_h,myGP_compact_cpp$varNoise)
    dist_mat[id_sensors,] <- dist_h
    
    write.table(hyp_mat, file = sprintf("exp_%s_capt_%s_hyperparameters.csv",           num2str(id_exp,fmt=0),num2str(list_sensors[id_sensors],fmt=0)),row.names = FALSE)
    write.table(dist_mat, file = sprintf("exp_%s_capt_%s_hyperparameter_distances.csv", num2str(id_exp,fmt=0),num2str(list_sensors[id_sensors],fmt=0)),row.names = FALSE)
    
    # now build another gp object such that you discard all the zero valued kernel functions, ie discard the informationless data points (Huygens' principle)
    coef(covWaveMaternCPP) <- coef(myGP_compact_cpp$covariance)
    var_h <- varVec(object = covWaveMaternCPP, X = wave_pos, compGrad = FALSE)
    
    id_kept_data <- which(var_h > 0)
    if (length(id_kept_data)==0){
      new_wave <- wave[1,] # will be zero in any case... just so it's not empty...   
    } else {
      new_wave <- wave[id_kept_data,]  
    }
    
    myGP_compact_cpp <- gp(formula = u ~ 1, data = new_wave, estim = FALSE,
                           cov = covWaveMaternCPP, varNoise = best_noise, parCovIni = parCovInitial, beta = 0)#, varNoise = 1e-5
    # print('Done!')
    print(sprintf('exp %s : running visual predictions...',num2str(id_exp,fmt=0)))
    predict_seq <- TRUE
    predict_for_errors <- TRUE
    if (predict_seq){ #compute a slice of the 3D solution for a visual estimation of the results
      print(sprintf('exp %s : running visual predictions...',num2str(id_exp,fmt=0)))
      ntot <- n_t*n_space^2
      pred_compact_mean_u <- cbind(myGrid_u,double(ntot)) # we use different prediction grids for u_0 and v_0
      pred_compact_std_u  <- cbind(myGrid_u,double(ntot))
      pred_compact_mean_v <- cbind(myGrid_v,double(ntot))
      pred_compact_std_v  <- cbind(myGrid_v,double(ntot))
      tic()
      par_current <- coef(myGP_compact_cpp$covariance)
      x0_u <- par_current[1:3]
      R_u  <- par_current[4]
      x0_v <- par_current[7:9]
      R_v  <- par_current[10]
      c    <- par_current[13]
      
      for (idt in 1:n_t){
        
        rows_h <- (1 + (idt-1)*n_space^2):(idt*n_space^2)
        grid_h_u <- myGrid_u[rows_h,]
        grid_h_v <- myGrid_v[rows_h,]
        
        # only predict values which are not known to be zero
        # first, select the points that need predictions
        print(sprintf('exp %s : visual predictions, u ; idt = %s...',num2str(id_exp,fmt=0),num2str(idt,fmt=0)))
        list_id_topredict_u <- selectPtsToPredict(grid_h_u,x0_u,R_u,c) # for recoding
        pts_topred_u <- grid_h_u[list_id_topredict_u,]
        
        if (nrow(pts_topred_u) == 0) {
          # then do nothing
        } else {
          pred_h_u <- predict(myGP_compact_cpp, newdata = pts_topred_u,seCompute = FALSE, covCompute = FALSE)
          rows_tofill_u <- list_id_topredict_u + (idt-1)*n_space^2
          
          pred_compact_mean_u[rows_tofill_u,5] <- pred_h_u$mean
          # pred_compact_std_u[rows_tofill,5]  <- pred_h_u$sd # uncomment to compute posterior standard deviation
        }
        
        print(sprintf('exp %s : visual predictions, v ; idt = %s...',num2str(id_exp,fmt=0),num2str(idt,fmt=0)))
        list_id_topredict_v <- selectPtsToPredict(grid_h_v,x0_v,R_v,c) # for recoding
        pts_topred_v <- grid_h_v[list_id_topredict_v,]
        if (nrow(pts_topred_v) == 0) {
          # then do nothing
        } else {
          pred_h_v <- predict(myGP_compact_cpp, newdata = pts_topred_v,seCompute = FALSE, covCompute = FALSE)
          rows_tofill_v <- list_id_topredict_v + (idt-1)*n_space^2
          
          pred_compact_mean_v[rows_tofill_v,5] <- pred_h_v$mean
          # pred_compact_std_v[rows_tofill,5]  <- pred_h_v$sd
        }
      }
      
      time_compact <- as.numeric(toc())
      print(sprintf('exp %s : saving visual predictions...',num2str(id_exp,fmt=0)))
      write.table(pred_compact_mean_u, file = sprintf("exp_%s_%ssensors_pred_wave_compact_Sim_u.csv",num2str(id_exp,fmt=0),num2str(nb_sensors_used,fmt=0)), row.names = FALSE)
      write.table(pred_compact_std_u,  file = sprintf("exp_%s_%ssensors_std_wave_compact_Sim_u.csv" ,num2str(id_exp,fmt=0),num2str(nb_sensors_used,fmt=0)), row.names = FALSE)
      write.table(pred_compact_mean_v, file = sprintf("exp_%s_%ssensors_pred_wave_compact_Sim_v.csv",num2str(id_exp,fmt=0),num2str(nb_sensors_used,fmt=0)), row.names = FALSE)
      write.table(pred_compact_std_v,  file = sprintf("exp_%s_%ssensors_std_wave_compact_Sim_v.csv" ,num2str(id_exp,fmt=0),num2str(nb_sensors_used,fmt=0)), row.names = FALSE)
    }
    
    if (predict_for_errors) {# compute kriging solution on 3D grid for 3D error computation
      print(sprintf('exp %s : running error predictions...',num2str(id_exp,fmt=0)))
      coord_min_u <- min(min(hyp_h[1:3] - hyp_h[4]),min(par_default[1:3]-par_default[4])) # make sure to compute the full space error and not just on some sub cube
      coord_max_u <- max(max(hyp_h[1:3] + hyp_h[4]),max(par_default[1:3]+par_default[4]))
      coord_min_v <- min(min(hyp_h[7:9] - hyp_h[10]),min(par_default[7:9]-par_default[10]))
      coord_max_v <- max(max(hyp_h[7:9] + hyp_h[10]),max(par_default[7:9]+par_default[10]))
      coord_min <- min(coord_min_u,coord_min_v)
      coord_max <- max(coord_max_u,coord_max_v)
      dx_error <- 1e-2 #Riemann sum space step
      n_error <- floor((coord_max-coord_min)/dx_error)
      n_err[id_sensors,1] <- n_error
      
      par_current <- coef(myGP_compact_cpp$covariance)
      x0_u <- par_current[1:3]
      R_u  <- par_current[4]
      x0_v <- par_current[7:9]
      R_v  <- par_current[10]
      c    <- par_current[13]
      
      tic()
      
      pred_compact_errors_mean <- c()
      
      for (idt in 1:n_t){
        
        ### Build a new regular grid for integral error computations
        grid_errors <- expand.grid(x = seq(coord_min, coord_max, length.out = n_error),
                                   y = seq(coord_min, coord_max, length.out = n_error),
                                   z = seq(coord_min, coord_max, length.out = n_error),
                                   t = dt*(idt-1), # t_rng = dt 
                                   KEEP.OUT.ATTRS = FALSE)
        
        pred_compact_errors_mean_h <- cbind(grid_errors,double(n_error^3))
        
        
        for (idz in 1:n_error){
          print(sprintf('Compact CPP for errors : idt = %s out of %s ; iteration %s out of %s',num2str(idt,fmt=0), num2str(n_t,fmt=0),num2str(idz,fmt=0), num2str(n_error,fmt=0)))
          rows_h <- (1 + (idz-1)*n_error^2):(idz*n_error^2)
          grid_h <- grid_errors[rows_h,]
          
          # only predict values which are not known to be zero
          # first, elect the points that need prediction
          
          list_id_topredict_u <- selectPtsToPredict(grid_h,x0_u,R_u,c) # for recoding
          pts_topred_u <- grid_h[list_id_topredict_u,]
          
          if (nrow(pts_topred_u) == 0) {
            # then do nothing
          } else {
            pred_h_u <- predict(myGP_compact_cpp, newdata = pts_topred_u,seCompute = FALSE, covCompute = FALSE)
            rows_tofill_u <- list_id_topredict_u + (idz-1)*n_error^2
            
            pred_compact_errors_mean_h[rows_tofill_u,5] <- pred_h_u$mean
          }
          
          list_id_topredict_v <- selectPtsToPredict(grid_h,x0_v,R_v,c) # for recoding
          pts_topred_v <- grid_h[list_id_topredict_v,]
          
          if (nrow(pts_topred_v) == 0) {
            # then do nothing
          } else {
            pred_h_v <- predict(myGP_compact_cpp, newdata = pts_topred_v,seCompute = FALSE, covCompute = FALSE)
            rows_tofill_v <- list_id_topredict_v + (idz-1)*n_error^2
            
            pred_compact_errors_mean_h[rows_tofill_v,5] <- pred_h_v$mean
          }
        }
        pred_compact_errors_mean <- rbind(pred_compact_errors_mean,pred_compact_errors_mean_h)
      }
      time_err[id_sensors] <- as.numeric(toc())
      print(sprintf('exp %s : saving predictions for error computations...',num2str(id_exp,fmt=0)))
      write.table(pred_compact_errors_mean, file = sprintf("exp_%s_%ssensors_pred_wave_compact_errors_Sim.csv",num2str(id_exp,fmt=0),num2str(nb_sensors_used,fmt=0)), row.names = FALSE)
      write.table(n_err, file = sprintf("exp_%s_n_error_vec.csv",num2str(id_exp,fmt=0)), row.names = FALSE)
    }
  }
  
  ### saving hyperparameters...
  print(sprintf('exp %s : saving all hyperparameters estimations...',num2str(id_exp,fmt=0)))
  write.table(hyp_mat, file = sprintf("exp_%s_hyperparameters.csv", num2str(id_exp,fmt=0)),row.names = FALSE)
  write.table(dist_mat, file = sprintf("exp_%s_hyperparameter_distances.csv", num2str(id_exp,fmt=0)),row.names = FALSE)
  write.table(n_err, file = sprintf("exp_%s_n_error_vec.csv",num2str(id_exp,fmt=0)), row.names = FALSE)
  write.table(time_hyp, file = sprintf("exp_%s_grp1_time_hyp.csv",num2str(id_exp,fmt=0)), row.names = FALSE)
  write.table(time_err, file = sprintf("exp_%s_grp1_time_err.csv",num2str(id_exp,fmt=0)), row.names = FALSE)
}
