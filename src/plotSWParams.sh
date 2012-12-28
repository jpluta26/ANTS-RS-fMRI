#!/usr/bin/env Rscript

# ---- about ----------
# john pluta 9/19/2012
# script to compute and statitistically compare small-world parameters in two groups of subjects
# INPUT: - two lists containing the explicit file path to the correlation matrices for all of the subjects included in each group.
# - number of time points in the fmri time series
# - specification of weather wavelet decomposition or standard correlation was used in generating the matrices
# OUTPUT: pdf file containing plots of the small world parameters (degree, clustering coefficient, path length, small worldness)
# and p-value of t-test between the two groups.
# ----------------------


# TODO: 
# -have option to print a pdf or a series of pngs


# ---- required libraries ----
library(methods)
library(brainwaver)
library(getopt) 
# ----------------------------














# -------------------------------------------------------------------------------------------- #
# ----------------------------------------- functions -----------------------------------------#
# -------------------------------------------------------------------------------------------- #

# ............ createAvgMat .......... #
# create matrix containing the average of all the correlation matrices for one group
# INPUT: list is a variable containing a list of paths to the correlation matrix for each subject
# in a group
# OUTPUT: a matrix containing the average of all correlation matrices
createAvgMat <- function(list)
{
	n.sub <- dim(list)[1];
	print(paste("Average correlation matrices of ", n.sub, " subjects... "))

	for(i in 1:n.sub)
	{
		# brainwaver adds header information to the correlation matrix
		# if using wavelet data, skip this
		cor.data <- read.table(list[i,]) 
		cor.data <- as.matrix(cor.data)
		
		
		
		# use the first subject to initialize some the storage matrix
		if(i == 1)
		{
			avgMat <- matrix(0, nrow=dim(cor.data)[1], ncol=dim(cor.data)[1])
		}
		
		# all subjects should have the same ROIs being tested, so correlation matrices should
		# have the same dimensions. 
		if( dim(avgMat)[1] != dim(cor.data)[1] )
		{
			print(paste("WARNING!"))
			print(paste(list[i,]))
			print(paste("Subject has dimensions that do not match the rest of the data.")) 
			print(paste("All subjects must have the same size correlation matrices. Exiting program."))
						 q(status=1);
		} else { avgMat <- avgMat + cor.data }
	}
	
	print(paste("Average matrix created"))
	return(avgMat / n.sub)
}
# ....................................... #



# .............. getSWParams ............... #
# function to compute small-world parameters from a given correlation matrix
# INPUT: correlation matrix (assumed to be group average), correlation threshold
# OUTPUT: dataframe containing small-world parameters
getSWParams <- function(avgMat, corthresh)
{
	# dataframe with small-world parameters
	params <- data.frame(in.degree.mean = 0, Cp.mean = 0, Lp.mean = 0, size.large.connex = 0, n.edges=0)
	n.nodes <- dim(avgMat)[1]
	
	# compute adjancency matrix - note that this includes postive and negative correlations
	adjMat <- const.adj.mat(avgMat, var.ind.mat=0, n.ind=0,
						thresh = 0.05, num.levels=1, sup=corthresh, 
						proc.length=n.tps, test.method="gaussian", use.tanh=FALSE)
	
	# since the adjacency matrix is symmetric, the number of edges is sum/2
	params$n.edges <- sum(adjMat)/2
	
	param.sw.brain <- small.world(adjMat, dat="reduced")
	in.degree.mean <- param.sw.brain$in.degree.mean
	in.Lp.mean <- param.sw.brain$Lp.mean
	in.Cp.mean <- param.sw.brain$Cp.mean
	
	degree.dist <- rowSums(adjMat)
		
	# create a random graph; this preserves degree distribution, number of nodes, and number of edges		
	# sink supresses the output		
	sink("/dev/null");
	rand <- equadist.rand.sw(10, dat="reduced", dist.mat=matrix(1,n.nodes,n.nodes), degree.dist=degree.dist);
	sink();
	
	# small-world parameters of the random graph
	rand.degree.mean <- rand$in.degree
	rand.Cp.mean <- rand$Cp.rand
	rand.Lp.mean <- rand$Lp.rand
	
	params$in.degree.mean <- in.degree.mean
	params$Lp.mean <- in.Lp.mean / rand.Lp.mean
	params$Cp.mean <- in.Cp.mean / rand.Cp.mean
	
	
	
	return(params)
}
# ...................................... #




# .......... plotparams ................ #
# function to plot different parameters
# INPUT:
plotparams <- function(ctlval, patval, y.label, main.label)
{
	print(paste("entering plotparams"))

	# determine y-range
	ymax <- max(ctlval, patval)
	ymax <- ymax + (ymax/10)  # add 10% extra y length to pad the graph
	cat("ymax = ", ymax)
	
	# plot patients versus controls
	plot(cor.seq, (1:n.cor)/2, type='n', xlab='Correlation threshold, R', 
	ylab=y.label, cex.axis=1, cex.lab=1, main=main.label, xlim=c(0,max(cor.seq)), ylim=c(0,ymax))
	

	# length of one increment
	inc <- cor.seq[2]-cor.seq[1]
	
	# setup plotting parameters
	plot.seq <- seq(from=0, to=(length(ctlval)-1)*inc, by=inc)

	lines(plot.seq, ctlval, type='l', col='blue', lwd=2)
	lines(plot.seq, patval, type='l', col='red', lwd=2)

	# do a t-test between the groups
	t <- t.test(ctlval, patval, alternative = "greater", paired=FALSE, var.equal=FALSE, conf.level=0.95)
	
	# attach the p-value to the plot
	mtext(sprintf(paste("p = %5.3f"), t$p.value));
	
}
# ...................................... #

# -------------------------------------------------------------------------------------------- #
# --------------------------------------- end functions -------------------------------------- #
# -------------------------------------------------------------------------------------------- #










# -------------------------------------------------------------------------------------------- #
# --------------------------------------------- user input ----------------------------------- #
# -------------------------------------------------------------------------------------------- #
Args <- commandArgs()
self <- Args[4]
self <- substring(self, 8, nchar(as.character(self)))


# define input options
spec = c(
	'help' , 'h', 0, "logical", "print help",
	'controllist', 'c', 1, "character", "file containing list of full path to invididual correlation matrices for control subjects",
	'patientlist', 'p', 1, "character", "file containing list of full paths to individual correlation matrices for patients",
    'controlfile', 'C', 1, "character", "file containing control data; use this to compare single subjects, or if you have already computed an average",
    'patientfile', 'P', 1, "character", "file containing patient data",
	'timepoints', 't', 1, "integer", "number of time-points the matrices were computed from",
	'figuretitle', 'f', 1, "character", "name to title the output figure (optional)",
	'outputfile',  'o', 1, "character", "name of the output file (optional)"
	)
	
# ADD HELP/USAGE!

spec = matrix(spec, ncol=5, byrow=TRUE)
opt = getopt(spec)

# .... help message
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
	ex <- paste(self," -t 120 -c ctrlist.txt -p patlist.txt -o swparams.pdf \n 
           self," -t 120 -C ctravgmat.txt -P patavgmat.txt \n\n")
	ex <- format(ex, width=length(ex), justify=c("left"))
	cat("\n")
	cat(ex)
	cat("Produces small world graph measures and stastical comparison of the two input groups, and writes output to pdf. \n")
	cat("The input files are text files specifying the full path to each individual subject's correlation matrix. \n")
	cat("The script makeSubList.sh can be used to generate these files.")
	q(status=1);
}

# ........... check for required input ............. #

if( is.null(opt$controlfile) & is.null(opt$controllist) )
{
    print(paste("Need to specify control list or control file"))
    print(paste("Exiting."))
    q(status=1)
}

if( is.null(opt$patientfile) & is.null(opt$patientlist))
{
    print(paste("Need to specify patient list or patient file"))
    print(paste("Exiting."))
    q(status=1)
}

if( is.null(opt$timepoints) )
{
	print(paste("Specify number of time points: -t "));
	print(paste("Exiting."))
	q(status=1);
}

# ................................................... #


# -------------------------------------------------------------------------------------------- #
# ------------------------------------------- end user input --------------------------------- #
# -------------------------------------------------------------------------------------------- #














# ============================================================================================= #
# ============================================================================================= #
# =============================================== MAIN ======================================== #
# ============================================================================================= #
# ============================================================================================= #
# setup constants

# names
defaultfilename="SW_params.pdf"
defaulttitle="Small World Paramteres"

# number of time-points that the correlation data was computed from
n.tps <- opt$timepoints


# sequence of correlations values to use as a threshold, to test parameters over a range of values
cor.seq <- seq(from=0, to=.4, by=.05)
n.cor <- length(cor.seq)

# empty vector to determine the size of the small-world dataframe
basevec <- rep(0, n.cor)
ctl.sm.params <- data.frame(in.degree.mean = basevec, Cp.mean = basevec, 
		Lp.mean = basevec, size.large.connex = basevec, n.edges = basevec)
pat.sm.params <- ctl.sm.params


if( !is.null(opt$controlfile) & !is.null(opt$patientfile))
{
    # read in the data lists
    ctrlist <- unname(as.matrix(read.table(opt$controllist)))
    patlist <- unname(as.matrix(read.table(opt$patientlist)))

    # create average matrices
    ctrAvgMat <- createAvgMat(ctrlist)
    patAvgMat <- createAvgMat(patlist)

    # might want to write these to .txt file?
    write.table(ctrAvgMat, file="ctrAvgMat.txt")
    write.table(patAvgMat, file="patAvgMat.txt")
} else
{
    
    ctrAvgMat <- as.matrix(read.table(eval(opt$controlfile)))
    patAvgMat <- as.matrix(read.table(eval(opt$patientfile)))
}

cat("read in the data")
# create a new index to have integer values for indexing
index <- 0

# populate small-world dataframes for controls and patients
for(cor.val in cor.seq)
{
	index <- index + 1
	
	ctl.sm.params[index,] <- getSWParams(ctrAvgMat, cor.val)
	
	print(paste("n.edges ", ctl.sm.params$n.edges[index]))

	# we expect patients to have a smaller number of edges for any given threshold
	# but we want to comparison to reflect the structure of the graphs, not just
	# be influenced by number of edges
	
	# this function choose a threshold to get the desired number of edges
	pat.cor.val <- choose.thresh.nbedges(patAvgMat, var.ind.mat = 0, n.ind = 0,
					   thresh = 0.05, nb.edges <- ctl.sm.params$n.edges[index],
					   test.method = "gaussian", proc.length = n.tps,
					   num.levels = 1, use.tanh = FALSE, max.iter = 50);
	
	
	pat.sm.params[index,] <- getSWParams(patAvgMat, pat.cor.val)
	

}





# ------------------------------------- plotting and statistical tests ------------------------------ #

if( is.null(opt$outputfile) )
{  outputfile=defaultfilename;  }  else  { outputfile = opt$outputfile }

if( is.null(opt$figuretitle) )
{  figuretitle=defaulttitle;  } else { figuretitle = opt$figuretitle }

# set up the pdf and plotting regions
pdf(outputfile, width=8.5, height=11)
par(mfrow=c(2,2))
par(mar=c(5, 4, 4, 2)) # inner margins
par(oma=c(1,1,3,1))    # outer margins

# vectors of the label names
y.label.vec <- c('k', parse(text="gamma"), parse(text="lambda"), parse(text="sigma") )
main.label.vec <- c("Mean Degree", "Mean Clustering Coefficient", "Mean Characteristic Path Length", "Small Worldness")


# plot all parameters
for(i in 1:4)
{
	print(paste("i = ", i))
	# for small worldness
	if( i == 4)
	{
		# vector of values
		ctlval <- ctl.sm.params[,2]/ctl.sm.params[,3]
		patval <- pat.sm.params[,2]/pat.sm.params[,3]
		
	} else 
	{  # for all others
		# vector of values
		ctlval <- ctl.sm.params[,i]
		patval <- pat.sm.params[,i]
				
		
	}
	
	# only retain valid values 
	ctlval <- ctlval[!(ctlval=="NaN")]
	ctlval <- ctlval[!(ctlval=="Inf")]
	patval <- patval[!(patval=="NaN")]
	patval <- patval[!(patval=="Inf")]

	smallest <- min(length(ctlval), length(patval))
	
	# need these to be the same length to plot against each other
    # shorten whichever vector is longer
	repeat {
		if( length(ctlval) > smallest )
		{
			ctlval <- ctlval[-(length(ctlval))]
		} else break;
	}

	repeat {
        if( length(patval) > smallest )
		{
			patval <- patval[-(length(patval))]
		} else break;
	}
		
	# setup labels
	y.label <- y.label.vec[i]
	main.label <- main.label.vec[i]
	
	# plot each parameter
	plotparams(ctlval, patval, y.label, main.label)
}



# figure caption
mtext(figuretitle, side=3, line=1, cex=1, col="black", outer=TRUE)  



dev.off()
# ---------------------------------- end plotting ---------------------------------------- #

# ============================================================================================= #
# ============================================================================================= #
# =========================================== END MAIN ======================================== #
# ============================================================================================= #
# ============================================================================================= #