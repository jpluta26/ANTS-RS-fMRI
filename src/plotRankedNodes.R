library(sna)
library(rgl)
library(brainwaver)
library(igraph)

# currently having major problems with X11 on mac
# try this on linux

n.tps = 240
corthresh=0.25
corMat <- read.table("/Users/pluta/Desktop/rsfmri_testdata/controls/n1/wave_cor_mat_level_1.txt", skip=3, header=FALSE)
corMat <- as.matrix(corMat)
adjMat <- const.adj.mat(corMat, var.ind.mat=0, n.ind=0, thresh = 0.05, num.levels=1, sup=corthresh, proc.length=n.tps, test.method="gaussian", use.tanh=FALSE)

g1 <- graph.adjacency(adjMat, mode=c("undirected") )
btwn <- betweenness(g1)

btwnRank <- rank(btwn)
nodeGrad <- colorRampPalette(c("red", "blue")) (10)
n.nodes = dim(corMat)[1]
nodeCol <- rep(" ", n.nodes)

for( i in 1:n.nodes) 
{
	if(btwnRank[i] >= 1 & btwnRank[i] <= 10)
	{ nodeCol[i] = nodeGrad[1] } else
if( btwnRank[i] >= 11 & btwnRank[i] <= 20)

{ nodeCol[i] = nodeGrad[2] } else
	if( btwnRank[i] >= 21 & btwnRank[i] <= 30)
	{ nodeCol[i] = nodeGrad[3] } else
	if( btwnRank[i] >= 31 & btwnRank[i] <= 40)
	{ nodeCol[i] = nodeGrad[4] } else
	if( btwnRank[i] >= 41 & btwnRank[i] <= 50)
	{ nodeCol[i] = nodeGrad[5] } else
	if( btwnRank[i] >= 51 & btwnRank[i] <=60 )
	{ nodeCol[i] = nodeGrad[6] } else
	if( btwnRank[i] >= 61 & btwnRank[i] <=70 )
	{ nodeCol[i] = nodeGrad[7] } else
	if( btwnRank[i] >= 71 & btwnRank[i] <=80 )
	{ nodeCol[i] = nodeGrad[8] } else
	if( btwnRank[i] >= 81 & btwnRank[i] <=90 )
	{ nodeCol[i] = nodeGrad[9] } else
	if( btwnRank[i] >= 91 & btwnRank[i] <=116 )
	{ nodeCol[i] = nodeGrad[10] } 
	
							
							 
}
coords=read.table("/Users/pluta/Desktop/aal_coords.txt")
coords <- as.matrix(coords)
gplot3d(adjMat, gmode="graph", diag=FALSE, vertex.col=nodeCol, new=FALSE, coord=coords)