library(igraph)
library(Matrix)


# Graphs ------------------------------------------------------------------
n <- 5000
tree <- make_tree(n,children=2,mode = c("undirected"))

scaleF <- barabasi.game(n,.05,directed=F)

erdosRenyi <- erdos.renyi.game(n,p=0.05,directed = F)

knMatrix <- matrix(data=1,n,n)
for(i in 1:n) knMatrix[i,i] <- 0
knMatrix

kn <- graph_from_adjacency_matrix(knMatrix,mode=c("undirected"))

graphs <- list(tree,scaleF,erdosRenyi,kn)
# Generate eigenvalues -----------------------------------------------------


#max(eigen(as_adjacency_matrix(erdosRenyi))$values)



g <- erdosRenyi
largestEigenvalue <- function(graph){
  M <- as_adj(graph)
  f2 <- function(x, extra=NULL) { cat("."); as.vector(M %*% x) }
  baev <- arpack(f2, sym=TRUE, options=list(n=vcount(graph), nev=1, ncv=20,
                                  which="LM", maxiter=200,ldv=0))
  max(abs(baev$values))
}

sapply(graphs,largestEigenvalue)
