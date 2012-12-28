#!/usr/bin/env Rscript

library(brainwaver)
library(igraph)

input is a correlation matrix

corMat <- read.table(file, header=FALSE, sep=" ")

# need to define n.tps and corthresh
adjMat <- const.adj.mat(corMat, var.ind.mat=0, n.ind=0,
						thresh = 0.05, num.levels=1, sup=corthresh, 
						proc.length=n.tps, test.method="gaussian", use.tanh=FALSE)
						
param.sw.brain <- small.world(adjMat, dat="reduced")

ewt <- c(0)
n.nodes = dim(corMat)[1]

for( x in c(1:n.nodes) )
{
	for( y in c(1:n.nodes) )
	{
		if( corMat[x ,y] > corthresh * corMat[x,y] < 1 ) { ewt<-c(ewt,im[x,y]) }
	}
}

ewt <- ewt[2:length(ewt)]

g1 <- graph.adjacency( adjMat, mode=c("undirected") )
g1 <- set.edge.attributes(g1, "weights", value=ewt*10)
btwn <- betweenness(g1)     # 1xN vec
ebtwn <- edge.betweenness(g1) # 1 x ewt vec
cl <- closeness(g1)   # 1xN vec
dg <- degree(g1)      # 1xN vec, can also get this from brainwaver for sure
hs <- hub.score(g1)$vector   # 1xN vec
pr <- page.rank(g1, weights=ewt)$vector  # 1xN vec
sp <- shortest.paths(g1, weights=1/ewt)  # NxN symmetric matrix
trans <- transitivity(g1, type="local")  # 1xN vec
						
names <- c("betweenness", "closeness", "degree", "hub score", "pagerank", "transitivity")

# format for output
params <- matrix(rbind(btwn, cl, dg, hs, pr, trans), nrow=6, ncol=n.nodes)
dimnames(params) <- list(names, seq(1:n.nodes) ) # this default, alterantively should be roi names