kalman <- function(A,B,H,R,Se,Phi,y){
  #Initialize the state vector at the stationary distribution
  T = nrow(y);
  l = ncol(y);
  n = nrow(Phi);
  s = matrix(0,T+1,n);
  P = array(0,c(T+1,n,n));
  s[1,] = matrix(0,1,n);
  a = solve(diag(n*n) - kronecker(Phi,Phi))%*%matrix(R%*%Se%*%Conj(t(R)),n*n,1);
  P[1,,] = matrix(a,n,n);
  # Kalman Filter Recursion
  sprime = matrix(0,n,1);
  Pprime = matrix(0,n,n);
  errorprediction = matrix(1,T,l);
  Varerrorprediction = array(1,c(T,l,l));
  liki = matrix(1,T,1);
  measurepredi = array(1,c(T,l));
  for (i in 1:T) {
    #updating step
    sprime = Phi%*%s[i,];
    Pprime = Phi%*%drop(P[i,,])%*%Conj(t(Phi)) + R%*%Se%*%Conj(t(R));
    # Prediction step
    yprediction = A + B%*%sprime;
    v = Conj(t(y[i,])) - yprediction;
    F = B%*%Pprime%*%Conj(t(B)) + H;
    kgain = Pprime%*%Conj(t(B))%*%solve(F);
    s[i+1,] = Conj(t(sprime + kgain%*%v));
    P[i+1,,] = Pprime - kgain%*%B%*%Pprime;
    errorprediction[i,] = Conj(t(v));
    Varerrorprediction[i,,] = F;
    liki[i] = -.5*l*log(2*pi) - .5*log(Det(F)) - .5*Conj(t(v))%*%solve(F)%*%v;
    measurepredi[i,] = as.matrix(y[i,] - t(v));
  }
  statepredi = s;
  varstatepredi = P;
  return(list(liki=liki,measurepredi=measurepredi,statepredi=statepredi,varstatepredi=varstatepredi))
}