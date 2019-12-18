prior <- function(Theta){
  a = matrix(0,1,4);
  b = matrix(1,1,4);
  P2 = 1/(b[1]-a[1]);
  P8 = 1/(b[2]-a[2]);
  P9 = 1/(b[3]-a[3]);
  P10 = 1/(b[4]-a[4]);
  # prior from gamma pdf
  para1 = c(2,1.5,0.5,0.5,7);
  para1 = matrix(para1,1,length(para1));
  para2 = c(0.5,0.25,0.25,0.5,2);
  para2 = matrix(para2,1,length(para2));
  b = para2^2/para1;
  a = para1/b;
  P1 = dgamma(Theta[1],shape=a[1],scale=b[1]);
  P3 = dgamma(Theta[3],shape=a[2],scale=b[2]);
  P4 = dgamma(Theta[4],shape=a[3],scale=b[3]);
  P5 = dgamma(Theta[5],shape=a[4],scale=b[4]);
  P6 = dgamma(Theta[6],shape = a[5],scale=b[5]);
  #Prior from normal pdf
  P7 = dnorm(Theta[7],.4,.2);
  # Prior from inverse gamma
  lnpdfig <- function(x,a,b){
    y = log(2) - log(gamma(b/2)) + (b/2)*log(b*a^2/2) - ((b+1)/2)*log(x^2) - b*a^2/(2*x^2);
    return(y)
  }
  P11 = exp(lnpdfig(Theta[11],0.4,4));
  P12 = exp(lnpdfig(Theta[12],1.00,4));
  P13 = exp(lnpdfig(Theta[13],0.50,4));
  f = c(P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13);
  f = prod(f);
  prior = log(f);
  return(prior)
}