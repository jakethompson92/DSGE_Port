dsgeliki <- function(para){
  T1 = model_solution(para)$T1;
  TC = model_solution(para)$TC;
  T0 = model_solution(para)$T0;
  TETA = model_solution(para)$TETA;
  RC = model_solution(para)$RC;
  retcode = model_solution(para)$retcode
  if(retcode==0){
    data = read.table('us.txt', header = FALSE, stringsAsFactors = FALSE);
    A = sysmat(T1,T0,para)$A;
    B = sysmat(T1,T0,para)$B;
    H = sysmat(T1,T0,para)$H;
    R = sysmat(T1,T0,para)$R;
    Se = sysmat(T1,T0,para)$Se;
    Phi = sysmat(T1,T0,para)$Phi;
    liki = kalman(A,B,H,R,Se,Phi,data)$liki;
    liki = sum(liki);
  }else{
    liki = -1E+12
  }
  return(liki)
}