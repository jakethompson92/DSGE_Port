objectiveunconstr <- function(Theta){
  prio = prior(Theta);
  if (prio==-Inf) {
    objective = -1e+12
  }else{
    liki = dsgeliki(Theta);
    objective = -(liki+prio);
  }
  return(objective)
}