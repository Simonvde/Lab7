#M is the adjacency matrix
#v is the vector of infected nodes (1->infected, 2->not infected)
#beta is the probability of infection, and gamma the probability of recovery

spread <- function(M,v,tmax,beta,gamma){
  n <- length(v)
  A <- matrix(nrow=tmax,ncol=n)
  A[1,] <- as.logical(v)
  for(i in 2:tmax){
    n_inf <- M%*%v  #a vector containing the number of infected neighbours of each node
    inf <- n_inf >= rgeom(n,beta)
    cured <- as.logical(rbinom(n=n,size=1,prob=gamma)*A[i-1,]) #vector infected of nodes which are cured (with probability gamma)
    A[i,] <- (A[i-1,] | inf) & !cured #an element is infected if it already was or gets infected this time step, and it hasn't been cured
  }
  return (A)
}