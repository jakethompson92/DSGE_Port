# Before running this, make sure to set the folder in which
# this script is placed as the present working directory
# This will enable R to source all and only the functions
# Written for this algorithm
sourceFolder <- function(folder, recursive = FALSE, ...) 
{ 
  files <- list.files(folder, pattern = "[.][rR]$", 
                      full.names = TRUE, recursive = recursive)
  if (!length(files))
    stop(simpleError(sprintf('No R files in folder "%s"', folder)))
  src <- invisible(lapply(files, source, ...))
  message(sprintf('%s files sourced from folder "%s"', length(src), folder))
}
sourceFolder('Sims')
library(tidyverse)
library(QZ)
library(complexplus)
library(MASS)
library(numDeriv)
# Declaring a parameter vector will make diagnostics much 
# easier
param = c(2.09,.6530,2.0,.65,.34,3.16,.51,1.5,4,2.5,.19,.65,.5);
param = matrix(param,1,length(param));
para = param
files = list.files()
files <- files[-c(2,3,4,6,7,8,13,14,16)]
lapply(files,source)
csminobj <- csminwelNew(objectiveconstr,t(para),diag(13),grad='numgrad',crit=1e-5,nit=200);
Theta = csminobj$xh
Theta[2] = exp(Theta[2])/(1+exp(Theta[2]));
Theta[8] = exp(Theta[8])/(1+exp(Theta[8]));
Theta[9] = exp(Theta[9])/(1+exp(Theta[9]));
Theta[10] = exp(Theta[10])/(1+exp(Theta[10]));
mode = t(Theta)
hess <- hessian(objectiveunconstr)
Sigma <- solve(hess) %>% Re() %>% matrix(,13,13)


#Metropolis Hastings algorithm

Nsim = 1000
c = .2;
c0 = .2
Nburn = as.integer(.5*Nsim)+2;
Thetasim = matrix(0,Nsim,13);
# Initialize by taking a draw from a distribution 
# centered at mode
go_on = 0;

while(go_on == 0){
  Thetac = mvrnorm(n=1,mu=mode,Sigma=c0*Sigma);
  go_on = (Thetac[8]<=1)*(Thetac[9]<=1)*(Thetac[10]<=1)*(Thetac[2]<=1);
}

Thetasim[1,] = Thetac

accept = 0;
obj <- dsgeliki(Thetasim[1,]) + prior(Thetasim[1,]);
obj <- Re(obj)
counter = 0;
logposterior = obj*matrix(1,Nsim,1);
#Sigma = proposal$Sigma



