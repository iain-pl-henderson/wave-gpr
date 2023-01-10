### clear workspace
rm(list = ls())

### indicate your workspace and library directory if required (e.g. if code running on server)
### you need to put the correct datasets in the same folder as this file
str_dir <- "/your_work_directory"
str_lib <- "/your_library_directory"

setwd(str_dir)

### uncomment to install the packages in the library directory
# install.packages("pracma",lib=str_lib)
# install.packages("lattice",lib=str_lib)
# install.packages("testthat",lib=str_lib)
# install.packages("nloptr",lib=str_lib)
# install.packages("kergp",lib=str_lib) # also requires RCPP

### import libraries
library("pracma",lib.loc=str_lib)
library("lattice",lib.loc=str_lib)
library("nloptr",lib.loc=str_lib)
library("testthat",lib.loc=str_lib)
library("kergp",lib.loc=str_lib)

### import dependencies
source('kernWaveMaternFun_pos.R') # contains WIGPR covariance function
source('selectPtsToPredict.R')    # small function which discards the informationless points

### Set booleans for the script
use_cpp <- 1

### new grid parameters 
n_exp <- 40
list_sensors <- c(3,5,10,15,20,25,30)
len_sensors <- length(list_sensors)

n_space <- 230
n_t <- 1
x_rng <- 1
t_rng <- 1e-7
dt <- t_rng/(n_t-1)

### Build manual kernels for the wave equation
print('Building manual kernels for the wave equation...')

par_default <- c(r_x = 0.5, r_y = 0.5, r_z = 0.5, R = 0.3, speed = 0.5, range = 0.03, var = 4) # good params for the cos ring

parLower_h <- c(r_x = 0, r_y = 0, r_z = 0, R = 0.03, speed = 0.2, range = 0.02, var = 0.1)
parUpper_h <- c(r_x = 1, r_y = 1, r_z = 1, R = 0.5, speed = 0.8, range = 2, var = 5)
parNames_h <- c("r_x", "r_y", "r_z", "R", "speed", "range", "var")

### Define the wave equation kernel for radially symmetric initial position
{
  inputNames_h <- c("x", "y", "z", "t")
  
  covWaveMaternCPP <- covMan(kernel = cppWave3D_pos,
                             acceptMatrix = FALSE,
                             hasGrad = FALSE,
                             label = "non stationary compactly supported wave kernel with cpp implementation",
                             d = 4,
                             inputs = inputNames_h,
                             parLower = parLower_h,
                             parUpper = parUpper_h,
                             parNames = parNames_h,
                             par = par_default) #c(speed = 1, sigma2 = 1, L = 3))
}

### Build a new regular grid for predictions
{
  myGrid <- expand.grid(x = seq(0, x_rng, length.out = n_space),
                        y = seq(0, x_rng, length.out = n_space),
                        z = 0.5, 
                        t = seq(0, t_rng, length.out = n_t), 
                        KEEP.OUT.ATTRS = FALSE)
}

for (id_exp in 1:n_exp){
  
  print(sprintf('exp %s out of %s',num2str(id_exp,fmt=0),num2str(n_exp,fmt=0)))
  
  ### Read data
  wave_full <- read.csv(sprintf("exp_%s_data_set_noisy.csv",num2str(id_exp,fmt=0)), sep = ",", header = FALSE)
  names(wave_full) <- c(inputNames_h, "u")
  nb_t_per_sensor <- 75
  
  ### Define data arrays for storing hyp result
  {
    hyp_mat <- matrix(0,len_sensors,8) # format : hyperparameters + estimated noise
    dist_mat <-matrix(0,len_sensors,3) # format : see hyperparameters
    n_err <- matrix(0,len_sensors,1) # list of integers needed for error computations
    time_hyp <- matrix(0,len_sensors)
    time_err <- matrix(0,len_sensors)
  }
  
  for (id_sensors in 1:len_sensors){
    
    nb_sensors_used <- list_sensors[id_sensors]
    
    print(sprintf('exp %s : nb sensors used : %s',num2str(id_exp,fmt=0),num2str(nb_sensors_used,fmt=0)))
    wave <- wave_full[1:(nb_sensors_used*nb_t_per_sensor),]
    wave_pos <- wave[,1:4]
    n_data <- nrow(wave_pos)
    
    ### use CPP source location kernel
    var_compact <- (0.45)^2 # put true noise value for sensibility study

    {      
      print(sprintf('exp %s : using exact hyperparameters...',num2str(id_exp,fmt=0)))
      
      parCovInitial <- par_default
      
      myGP_compact_cpp <- gp(formula = u ~ 1, data = wave, estim = FALSE, varNoise = var_compact,
                             cov = covWaveMaternCPP, parCovIni = parCovInitial, beta = 0)
    }
    
    hyp_h <- coef(myGP_compact_cpp$covariance)
    # for cos
    dist_x0 <- sqrt(sum((hyp_h[1:3] - par_default[1:3])^2))
    dist_R <- abs(hyp_h[4]-par_default[4])
    dist_c <- abs(hyp_h[5]-par_default[5])

    dist_h <-c(dist_x0,dist_R,dist_c)
    
    # Store estimated hyperparameters (useless for sensibility analysis)
    hyp_mat[id_sensors,] <- c(hyp_h,myGP_compact_cpp$varNoise)
    dist_mat[id_sensors,] <- dist_h
    
    # now build another gp object such that you discard all the zero valued kernel functions, ie discard the informationless data points
    coef(covWaveMaternCPP) <- coef(myGP_compact_cpp$covariance)
    var_h <- varVec(object = covWaveMaternCPP, X = wave_pos, compGrad = FALSE)
    
    id_kept_data <- which(var_h > 0)
    if (length(id_kept_data)==0){
      new_wave <- wave[1,] # will be zero in any case... just so it's not empty...   
    } else {
      new_wave <- wave[id_kept_data,]  
    }
    
    myGP_compact_cpp <- gp(formula = u ~ 1, data = new_wave, estim = FALSE,
                           cov = covWaveMaternCPP, varNoise = var_compact, parCovIni = parCovInitial, beta = 0)#, varNoise = 1e-5
    # print('Done!')
    print(sprintf('exp %s : running visual predictions...',num2str(id_exp,fmt=0)))
    predict_seq <- TRUE # parallel computation not implemented
    predict_for_errors <- TRUE
    if (predict_seq){
      
      ntot <- n_t*n_space^2
      pred_compact_mean <- cbind(myGrid,double(ntot))
      pred_compact_std <- cbind(myGrid,double(ntot))
      tic()
      par_current <- coef(myGP_compact_cpp$covariance)
      x0 <- par_current[1:3]
      R <- par_current[4]
      c <- par_current[5]
      
      for (idt in 1:n_t){
        
        print(sprintf('Compact CPP : iteration %s out of %s',num2str(idt,fmt=0), num2str(n_t,fmt=0)))
        rows_h <- (1 + (idt-1)*n_space^2):(idt*n_space^2)
        grid_h <- myGrid[rows_h,]
        
        # only predict values which are not known to be zero
        # first, select the points that need prediction
        list_id_topredict <- selectPtsToPredict(grid_h,x0,R,c) # for recoding
        pts_topred <- grid_h[list_id_topredict,]
        
        if (nrow(pts_topred) == 0) {
          # then do nothing
        } else {
          pred_h <- predict(myGP_compact_cpp, newdata = pts_topred)
          rows_tofill <- list_id_topredict + (idt-1)*n_space^2
          
          pred_compact_mean[rows_tofill,5] <- pred_h$mean
          pred_compact_std[rows_tofill,5] <- pred_h$sd
        }
      }
      
      time_compact <- toc()
      print(sprintf('exp %s : saving visual predictions...',num2str(id_exp,fmt=0)))
      write.table(pred_compact_mean, file = sprintf("exp_%s_%ssensors_pred_wave_compact_Sim.csv",num2str(id_exp,fmt=0),num2str(nb_sensors_used,fmt=0)), row.names = FALSE)
      write.table(pred_compact_std, file = sprintf("exp_%s_%ssensors_std_wave_compact_Sim.csv",num2str(id_exp,fmt=0),num2str(nb_sensors_used,fmt=0)), row.names = FALSE)
      print('Done!')
    }
    
    if (predict_for_errors) { # compute kriging solution on 3D grid for 3D error computation
      coord_min <- min(min(hyp_h[1:3] - hyp_h[4]),min(par_default[1:3]-par_default[4])) # make sure to compute the full space error and not just on some sub cube
      coord_max <- max(max(hyp_h[1:3] + hyp_h[4]),max(par_default[1:3]+par_default[4]))
      dx_error <- 1e-2
      n_error <- floor((coord_max-coord_min)/dx_error)
      n_err[id_sensors,1] <- n_error
      
      par_current <- coef(myGP_compact_cpp$covariance)
      x0 <- par_current[1:3]
      R <- par_current[4]
      c <- par_current[5]
      tic()
      
      ### Build a new regular grid for integral error computations
      grid_errors <- expand.grid(x = seq(coord_min, coord_max, length.out = n_error),
                                 y = seq(coord_min, coord_max, length.out = n_error),
                                 z = seq(coord_min, coord_max, length.out = n_error),
                                 t = t_rng, # t_rng = dt 
                                 KEEP.OUT.ATTRS = FALSE)
      
      pred_compact_errors_mean <- cbind(grid_errors,double(n_error^3))
      
      for (idz in 1:n_error){
        print(sprintf('Compact CPP for errors : iteration %s out of %s',num2str(idz,fmt=0), num2str(n_error,fmt=0)))
        rows_h <- (1 + (idz-1)*n_error^2):(idz*n_error^2)
        grid_h <- grid_errors[rows_h,]
        #grid_h <- slice(grid_errors,rows_h)
        
        # only predict values which are not known to be zero
        # first, elect the points that need prediction
        list_id_topredict <- selectPtsToPredict(grid_h,x0,R,c) # for recoding
        pts_topred <- grid_h[list_id_topredict,]
        
        if (nrow(pts_topred) == 0) {
          # then do nothing
        } else {
          pred_h <- predict(myGP_compact_cpp, newdata = pts_topred)
          rows_tofill <- list_id_topredict + (idz-1)*n_error^2
          
          pred_compact_errors_mean[rows_tofill,5] <- pred_h$mean
        }
      }
      
      time_err[id_sensors] <- toc()
      print(sprintf('exp %s : saving predictions for error computations...',num2str(id_exp,fmt=0)))
      write.table(pred_compact_errors_mean, file = sprintf("exp_%s_%ssensors_pred_wave_compact_errors_Sim.csv",num2str(id_exp,fmt=0),num2str(nb_sensors_used,fmt=0)), row.names = FALSE)
      print('Done!')
      
    }
  }
  
  ### saving hyperparameters...
  print(sprintf('exp %s : saving all hyperparameters estimations...',num2str(id_exp,fmt=0)))
  write.table(hyp_mat, file = sprintf("exp_%s_hyperparameters.csv", num2str(id_exp,fmt=0)),row.names = FALSE)
  write.table(dist_mat, file = sprintf("exp_%s_hyperparameter_distances.csv", num2str(id_exp,fmt=0)),row.names = FALSE)
  write.table(n_err, file = sprintf("exp_%s_n_error_vec.csv",num2str(id_exp,fmt=0)), row.names = FALSE)
  write.table(time_hyp, file = sprintf("exp_%s_time_hyp.csv",num2str(id_exp,fmt=0)), row.names = FALSE)
  write.table(time_err, file = sprintf("exp_%s_time_err.csv",num2str(id_exp,fmt=0)), row.names = FALSE)
}
