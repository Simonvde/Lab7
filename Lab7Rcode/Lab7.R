library(igraph)
library(Matrix)


# Infection Function ------------------------------------------------------

#M is the adjacency matrix
#v is the vector of infected nodes (1->infected, 0->not infected)
#beta is the probability of infection, and gamma the probability of recovery

spread <- function(M,v,tmax,beta,gamma){
  n <- length(v)
  A <- matrix(nrow=tmax,ncol=n)
  A[1,] <- as.logical(v)
  for(i in 2:tmax){
    n_inf <- M%*%v  #a vector containing the number of infected neighbours of each node
    inf <- n_inf > rgeom(n,beta)
    cured <- as.logical(rbinom(n=n,size=1,prob=gamma)*A[i-1,]) #vector infected of nodes which are cured (with probability gamma)
    A[i,] <- as.logical((A[i-1,] | inf) & !cured) #an element is infected if it already was or gets infected this time step, and it hasn't been cured
  }
  return (A)
}

# Graphs ------------------------------------------------------------------
n <- 2000
tree <- make_tree(n,children=2,mode = c("undirected"))

scaleF <- barabasi.game(n,.05,directed=F)

erdosRenyi <- erdos.renyi.game(n,p=0.05,directed = F)

knMatrix <- matrix(data=1,n,n)
for(i in 1:n) knMatrix[i,i] <- 0


kn <- graph_from_adjacency_matrix(knMatrix,mode=c("undirected"))

graphs <- list(tree,scaleF,erdosRenyi,kn)
# Generate eigenvalues -----------------------------------------------------


#max(eigen(as_adjacency_matrix(erdosRenyi))$values)
largestEigenvalue <- function(graph){
  M <- as_adj(graph)
  f2 <- function(x, extra=NULL) { cat("."); as.vector(M %*% x) }
  baev <- arpack(f2, sym=TRUE, options=list(n=vcount(graph), nev=1, ncv=20,
                                  which="LM", maxiter=200,ldv=0))
  max(abs(baev$values))
}

sapply(graphs,largestEigenvalue)



# Simulate -------------------------------------------------------------------


infected <- sample(c(1,0),n,replace=T,prob=c(0.05,0.95))
as.logical(infected)
treeSpread <- spread(as_adj(kn),v=infected,tmax=100,beta=.1,gamma=.05)

infectedEvolution <- function(spread){
  infected <- c()
  for(i in 1:length(spread[,1])){
    infected[i] <- length(which(spread[i,]==T))
  }
  infected
}
infectedEvolution(treeSpread)

plot(infectedEvolution(treeSpread))

simulate <- function(graph,beta=.4,gamma=.8,p0=.05,tmax=30){
  initialInfected <- sample(c(1,0),length(V(graph)),replace=T,prob=c(p0,1-p0))
  treeSpread <- spread(as_adj(graph),v=initialInfected,tmax=tmax,beta=beta,gamma=gamma)
  infectedEvolution(treeSpread)
}


# beta,gamma close to threshold -------------------------------------------

epsilon=.05
thresholds = 1/sapply(graphs,largestEigenvalue)
gamma <- .4
beta <- gamma*thresholds[1]+epsilon
yvalues <- simulate(tree,p0=.05,beta=beta,gamma=gamma,tmax=30)
plot(yvalues,ylim = c(0,max(yvalues)))
beta <- gamma*thresholds[1]-epsilon
lines(simulate(tree,p0=.05,beta=beta,gamma=gamma,tmax=30))
lines((beta-gamma)*n/beta)
