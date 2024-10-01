# wave_gpr
This code can be used to reproduce the numerical results of the paper "Covariance models and Gaussian process regression for the wave equation. Application to related inverse problems" (I. Henderson, P. Noble, O. Roustant), In: Journal of Computational Physics 494(2023), Paper No. 112519,  available at https://arxiv.org/abs/2311.05205 .

Specifically, it contains
1) a MATLAB FDTD (Finite Difference Time Domain) script for solving the 3D free-field wave equation, located in wave_equation_simulator/FDTD_wave_eq_3D_git.m .
2) several R scripts corresponding to each numerical experiment, located in wave_equation_kriging/code .

For example, for performing Kriging on one solution of the wave equation with non zero initial position and velocity, the R script for launching the numerical experiment is located in wave-gpr/wave_equation_kriging/code/mix/multistart/GP_wave_3D_multistart_mix.R .
This script performs Kriging on a given numerical solution to the wave equation, using the R package kergp. The corresponding wave equation-tailored covariance function plugged in kergp is coded in kernWaveMaternFun_pos_spd.R, using a RCpp implementation.

In the case of the analysis of the sensibility of our method w.r.t. the sensor locations, the corresponding R script is located in wave_equation_kriging/code/mix/sensibility/GP_wave_3D_sensib_mix.R .

More details on each experiment are given in the article.

