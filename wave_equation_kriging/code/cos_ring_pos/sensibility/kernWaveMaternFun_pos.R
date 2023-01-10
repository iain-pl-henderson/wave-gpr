
cppFunction('

  NumericVector cppWave3D_pos(NumericVector x1, NumericVector x2,
                                      NumericVector par);

  NumericVector phitan(double x);
                                 
  NumericVector cppKernMatern5_2(double x1, double x2,
                                 NumericVector par);
            
  NumericVector cppKernMatern5_2(double x1, double x2,
                                 NumericVector par){
    NumericVector K(1);
    double h, z, emz;
    h = (x1 - x2) / par[0];
    z = sqrt(5) * fabs(h);
    emz = exp(-z);
    K[0] = par[1] * (1.0 + z * ( 1.0 + z / 3.0) ) * emz;
    return(K);
  }                                    
  
  NumericVector phitan(double r){
    double pi = 3.14159265359;
    double alpha = 0.95;
    
    NumericVector val(1);
    if (r < alpha){
      val[0] = 1;
    } else if(r >= 1){
      val[0] = 0;
    } else {
      double th = std::tan(pi*0.5*(r-alpha)/(1-alpha));
      val[0] = exp(-1/(th*th));
    }
    return val;
  }
  
  
  
  NumericVector cppWave3D_pos(NumericVector x1, NumericVector x2,
                                      NumericVector par){
                                      
  NumericVector x1_loc(x1.begin(), x1.end());
  NumericVector x2_loc(x2.begin(), x2.end());
  
  double x1n, x2n, s1x1, s2x1, s1x2, s2x2, R;
  NumericVector k11(1), k12(1), k21(1), k22(1);
  NumericVector K(1), parK0(2);
  R = par[3];
  
  for (int i=0; i<3; i++) {
    x1_loc[i] = x1_loc[i] - par[i];
    x2_loc[i] = x2_loc[i] - par[i];
  }
  
  x1_loc[3] = x1_loc[3] * par[4];
  x2_loc[3] = x2_loc[3] * par[4];
  
  x1n = sqrt(x1_loc[0] * x1_loc[0] + x1_loc[1] * x1_loc[1] + x1_loc[2] * x1_loc[2]);

  s1x1 = x1n - x1_loc[3];
  
  if (fabs(s1x1) > R) {
    K[0] = 0;
    return(K);
  }
  
  x2n = sqrt(x2_loc[0] * x2_loc[0] + x2_loc[1] * x2_loc[1] + x2_loc[2] * x2_loc[2]);
  
  s1x2 = x2n - x2_loc[3];
  
  if (fabs(s1x2) > R) {
    K[0] = 0;
    return(K);
  }
  
  s2x1 = x1n + x1_loc[3];
  s2x2 = x2n + x2_loc[3];
  
  parK0[0] = par[5];
  parK0[1] = par[6];
  
  if (s2x1 > R){
  
    k22 = 0;
    k21 = 0;
    
    if (s2x2 > R){
      k12 = 0;
    } else {
      k12 = s1x1 * s2x2 * cppKernMatern5_2(s1x1*s1x1, s2x2*s2x2, parK0) * phitan(s1x1/R) * phitan(s2x2/R);
    }
    
  } else {
  
      if (s2x2 > R){
        k22 = 0;
        k12 = 0;
      } else {
        k12 = s1x1 * s2x2 * cppKernMatern5_2(s1x1*s1x1, s2x2*s2x2, parK0)  * phitan(s1x1/R) * phitan(s2x2/R);
        k22 = s2x1 * s2x2 * cppKernMatern5_2(s2x1*s2x1, s2x2*s2x2, parK0)  * phitan(s2x1/R) * phitan(s2x2/R);
      }
    k21 = s2x1 * s1x2 * cppKernMatern5_2(s2x1*s2x1, s1x2*s1x2, parK0)  * phitan(s2x1/R) * phitan(s1x2/R);
  }
  
  k11 = s1x1 * s1x2 * cppKernMatern5_2(s1x1*s1x1, s1x2*s1x2, parK0)  * phitan(s1x1/R) * phitan(s1x2/R);
  
  K = (k11 + k22 + k12 + k21)/ (4 * x1n * x2n );
  
  return(K);
  }')

# if (s2x2 > R){
#   k22 = 0;
#   k12 = 0;
# } else {
#   k12 = s1x1 * s2x2 * cppKernMatern5_2(s1x1*s1x1, s2x2*s2x2, parK0);
# }