sysmat <- function(T1,T0,para){
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
  eq_y = 1;
  eq_pi = 2;
  eq_ffr = 3;
  ny = 3;
  y_t = 1;
  pi_t = 2;
  R_t = 3;
  y1_t = 4;
  g_t = 5;
  z_t = 6;
  Ey_t1 = 7;
  Epi_t1 = 8;
  z_sh = 1;
  g_sh = 2;
  R_sh = 3;
  nep = ncol(T0);
  Phi = T1;
  R = T0;
  Se = matrix(0,nep,nep);
  Se[z_sh,z_sh] = sigma_z^2;
  Se[g_sh,g_sh] = (sigma_g)^2;
  Se[R_sh,R_sh] = (sigma_R)^2;
  A = matrix(0,ny,1);
  A[eq_y,1] = gammaQ;
  A[eq_pi,1] = piA;
  A[eq_ffr,1] = piA + rA + 4*gammaQ;
  nstate = ncol(Phi);
  B = matrix(0,ny,nstate);
  B[eq_y,y_t] = 1;
  B[eq_y,y1_t] = -1;
  B[eq_y,z_t] = 1;
  B[eq_pi,pi_t] = 4;
  B[eq_ffr,R_t] = 4;
  H = matrix(0,ny,ny);
  H[eq_y,y_t] = (0.20*0.579923)^2;
  H[eq_pi,pi_t] = (0.20*1.470832)^2;
  H[eq_ffr,R_t] = (0.20*2.237937)^2;
  return(list(A=A,B=B,H=H,R=R,Se=Se,Phi=Phi))
}