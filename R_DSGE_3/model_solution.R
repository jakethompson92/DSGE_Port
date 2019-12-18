model_solution <- function(para){
  tau = para[1];
  kappa = para[2];
  psi1 = para[3];
  psi2 = para[4];
  rA = para[5];
  piA = para[6];
  gammaQ = para[7];
  rho_R = para[8];
  rho_g = para[9];
  rho_z = para[10];
  sigma_R = para[11];
  sigma_g = para[12];
  sigma_z = para[13];
  bet = 1/(1+rA/400);
  # Solve DSGE model
  retcode = 0;
  valid = 1;
  # equation indices
  eq_1 = 1;
  eq_2 = 2;
  eq_3 = 3;
  eq_4 = 4;
  eq_5 = 5;
  eq_6 = 6;
  eq_7 = 7;
  eq_8 = 8;
  # variable indices
  y_t = 1;
  pi_t = 2;
  R_t = 3;
  y1_t = 4;
  g_t = 5;
  z_t = 6;
  Ey_t1 = 7;
  Epi_t1 = 8;
  # Expectation error indices
  ey_sh = 1;
  epi_sh = 2;
  # Shock indices
  z_sh = 1;
  g_sh = 2;
  R_sh = 3;
  # Summary 
  neq = 8;
  neta = 2;
  neps = 3;
  # Initialize your matrices 
  GAM0 = matrix(0,neq,neq);
  GAM1 = matrix(0,neq,neq);
  C = matrix(0,neq,1);
  PSI = matrix(0,neq,neps);
  PPI = matrix(0,neq,neta);
  # Equilibrium conditions: Canonical system
  GAM0[eq_1,y_t] = 1;
  GAM0[eq_1,R_t] = 1/tau;
  GAM0[eq_1,g_t] = -(1-rho_g);
  GAM0[eq_1,z_t] = -rho_z/tau;
  GAM0[eq_1,Ey_t1] = -1;
  GAM0[eq_1,Epi_t1] = -1/tau;
  # 2 
  GAM0[eq_2,y_t] = -kappa;
  GAM0[eq_2,pi_t] = 1;
  GAM0[eq_2,g_t] = kappa;
  GAM0[eq_2,Epi_t1] = -bet;
  #3 
  GAM0[eq_3,y_t] = -(1-rho_R)*psi2;
  GAM0[eq_3,pi_t] = -(1-rho_R)*psi1;
  GAM0[eq_3,R_t] = 1;
  GAM0[eq_3,g_t] = (1-rho_R)*psi2;
  GAM1[eq_3,R_t] = rho_R;
  PSI[eq_3,R_sh] = 1;
  #4 
  GAM0[eq_4,y1_t] = 1;
  GAM1[eq_4,y_t] = 1;
  #5 
  GAM0[eq_5,g_t] = 1;
  GAM1[eq_5,g_t] = rho_g;
  PSI[eq_5,g_sh] = 1;
  #6 
  GAM0[eq_6,z_t] = 1;
  GAM1[eq_6,z_t] = rho_z;
  PSI[eq_6,z_sh] = 1;
  #7
  GAM0[eq_7,y_t] = 1;
  GAM1[eq_7,Ey_t1] = 1;
  PPI[eq_7,ey_sh] = 1;
  #8 
  GAM0[eq_8,pi_t] = 1;
  GAM1[eq_8,Epi_t1] = 1;
  PPI[eq_8,epi_sh] = 1;
  # QZ(Generalized Schur) decomposition by GENSYS
  T1 = gensys(GAM0,GAM1,C,PSI,PPI,1+1e-8)$G1;
  TC = gensys(GAM0,GAM1,C,PSI,PPI,1+1e-8)$C;
  T0 = gensys(GAM0,GAM1,C,PSI,PPI,1+1e-8)$impact;
  TY = gensys(GAM0,GAM1,C,PSI,PPI,1+1e-8)$fmat;
  M = gensys(GAM0,GAM1,C,PSI,PPI,1+1e-8)$fwt;
  TZ = gensys(GAM0,GAM1,C,PSI,PPI,1+1e-8)$ywt;
  TETA = gensys(GAM0,GAM1,C,PSI,PPI,1+1e-8)$gev;
  GEV = gensys(GAM0,GAM1,C,PSI,PPI,1+1e-8)$eu;
  RC = gensys(GAM0,GAM1,C,PSI,PPI,1+1e-8)$loose;
  return(list(T1=T1,TC=TC,T0=T0,TETA=TETA,RC=RC,retcode=retcode))
}