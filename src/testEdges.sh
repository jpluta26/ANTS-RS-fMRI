#!/usr/bin/env Rscript

## edit this
# ---- about ----------
# john pluta 9/19/2012
# script to compute and statitistically compare small-world parameters in two groups of subjects
# INPUT: - two lists containing the explicit file path to the correlation matrices for all of the subjects included in each group.
# - number of time points in the fmri time series
# - specification of weather wavelet decomposition or standard correlation was used in generating the matrices
# OUTPUT: pdf file containing plots of the small world parameters (degree, clustering coefficient, path length, small worldness)
# and p-value of t-test between the two groups.
# ----------------------



# ---- required libraries ----
library(methods)
library(brainwaver)
library(getopt) 
# ----------------------------



createCorMat3D <- function(list)
{
	n.sub <- dim(list)[1];
	
	# read in the first subject to initialize some variables
	if( ! is.null(opt$wavelet) )
		 { cor.data <- read.table(list[1,], skip=3) }
		else  {  cor.data <- read.table(list[1,]) }

	n.nodes <- dim(cor.data)[1]
	
	# create a 3D matrix to store the distribution of correlation values for each node
	dataMat <- array(0, dim=c(n.nodes, n.nodes, n.sub) )
	dataMat[,,1] <- cor.data
	
	for(i in 2:n.sub)
	{
		# brainwaver adds header information to the correlation matrix
		# if using wavelet data, skip this
		if( ! is.null(opt$wavelet) )
		 { cor.data <- read.table(list[i,], skip=3) }
		else  {  data <- read.table(list[i,]) }

		dataMat[,,i] <- as.matrix(cor.data)
	}
}
		
		
		
		
oneSidePvalue <- function(dataMat)
{	
	# get number of nodes
	n.nodes=dim(dataMat)[1]
	
	# create an empty matrix to store p-values of the one-sided t-tests
	oneside.pval <- matrix(1, nrow=n.nodes, ncol=n.nodes)
	
	
	# the correlation matrix is symmetric and the diagonals are 1. so we only need to test off-diagonal
	# and upper-traingular values. first create an index for these:
	index <- upper.tri(dataMat[,,1]) # logical vector containing TRUE for all upper triangular positions
	
	# loop through the index matrix
	for(i in 1:n.nodes)
	{
		for(j in 1:n.nodes)
		{
			# for valid indices, retrive the entire distribution of values and test versus mu=0
			# store output in p-value matrix
			if( index[i,j] == TRUE)
			{
				x <- dataMat[i,j,];
				t <- t.test(x, alternative=c("two.sided"), mu=0)
				oneside.pval[i,j] <- t$p.value
		}
	}
	
	# return the matrix of p-values
	return(oneside.pval)
	
}		


# --------------------------------- user input ----------------------------------- #
Args <- commandArgs()
self <- Args[4]
self <- substring(self, 8, nchar(as.character(self)))


# define input options
spec = c(
	'help' , 'h', 0, "logical", "print help",
	'controlfile', 'c', 1, "character", "file containing list of full path to invididual correlation matrices for control subjects",
	'patientfile', 'p', 1, "character", "file containing list of full paths to individual correlation matrices for patients",
	'timepoints', 't', 1, "integer", "number of time-points the matrices were computed from",
	'wavelet', 'w', 0, "logical", "matrices were created by wavelet decomposition",
	'correlation', 'r', 0, "logical", "matrices were created by correlation",
	'figuretitle', 'f', 1, "character", "name to title the output figure (optional)",
	'outputfile',  'o', 1, "character", "name of the output file (optional)"
	)
	
# ADD HELP/USAGE!

spec = matrix(spec, ncol=5, byrow=TRUE)
opt = getopt(spec)

# .... help message ## EDIT THIS
if ( !is.null(opt$help) || length(opt) == 1) 
{
	cat("\n")
	cat(paste(self,"\n"))
	
	for( x in 1:nrow(spec) )
	{
		cat("\n")
		longopt <- paste("--", spec[x,1], sep='')
		shortopt <- paste("-", spec[x,2], sep='')
		hlist <- paste(shortopt, "|", longopt, spec[x,5], "\n \n")
		cat(format(hlist, width=40, justify = c("left")))
	}
	
	cat(format("Example: \n", width=40, justify = c("left")))
	ex <- paste(self," -w -t 120 -c ctrlist.txt -p patlist.txt -o swparams.pdf \n \n")
	ex <- format(ex, width=length(ex), justify=c("left"))
	cat("\n")
	cat(ex)
	cat("Produces small world graph measures and stastical comparison of the two input groups, and writes output to pdf. \n")
	cat("The input files are text files specifying the full path to each individual subject's correlation matrix. \n")
	cat("The script makeSubList.sh can be used to generate these files.")
	q(status=1);
}

# ........... check for required input ............. #
if( is.null(opt$controlfile) )
{
	print(paste("Specify file containing path to control matrices: -c"))
	print(paste("Exiting."))
	q(status=1);
}

if( is.null(opt$patientfile) )
{
	print(paste("Specify file containing path to control matrices: -p"))
	print(paste("Exiting."))
	q(status=1);
}

if( is.null(opt$timepoints) )
{
	print(paste("Specify number of time points: -t "));
	print(paste("Exiting."))
	q(status=1);
}

if( is.null(opt$wavelet) & is.null(opt$correlation) )
{
	print(paste("Identify matrices as either wavelet decomposition or standard correlation (-w or -r)"))
	print(paste("Exiting."))
	q(status=1);
} else
if ( !(is.null(opt$wavelet)) & !(is.null(opt$correlation)) )
{
	print(paste("Identify matrices as either wavelet decomposition or standard correlation (-w or -r)"))
	print(paste("Exiting."))
	q(status=1);
} 
# ................................................... #



# ------------------------------------------- end user input --------------------------------- #







# =============================================== MAIN ======================================== #
# setup constants
n.tps <- opt$timepoints


# read in the data lists
ctrlist <- unname(as.matrix(read.table(opt$controlfile)))
ctrDataMat <- createCorMat3D(ctrlist)

patlist <- unname(as.matrix(read.table(opt$patientfile)))
patDataMat <- createCorMat3D(patlist)

# test each edge versus 0
ctl.oneside.pval <- oneSidePval(ctrDataMat)
pat.oneside.pval <- oneSidePval(patDataMat)

# one sample t-test vs 0 in EITHER group must be signficant to show that the correlation is not spurious
# then test the two samples against each other
n.nodes <- dim(ctl.oneside.pval)[1]

# matrix to store the p-values of the control vs patient tests
grp.test.pval <- matrix(1, nrow=n.nodes, ncol=n.nodes)


for(i in 1:n.nodes)
{
	for(j in 1:n.nodes)
	{	
		# if either is significant, test against each other
		if( (ctl.oneside.pval[i,j] <= 0.05) | (pat.oneside.pval[i,j] <= 0.05) )
		{
			x <- ctrDataMat[i,j,]
			y <- patDataMat[i,j,]
			t <- t.test(x, y, alternative=c("two.sided"), paired=FALSE, var.equal=FALSE, conf.level=0.95)
			grp.test.pval[i,j] <- t$p.value
		}
	}
}

# how to visualize the results?
# could do a 3D plot
# textfile summary
