#!/usr/bin/env Rscript

# --------- about --------
# avants and pluta 2012
# some code is taken from the brainwaver manual, by sophie achard
#
# script to perform filtering of time-series data and create correlation matrix/matrices
# uses either standard correlation within a specified filtering range, or wavelet decomposition with 6 levels
# INPUT: fmri time-series, nuisance parameters, motion-correction parameters (.csv files)
# number of dummy scans, tr, and optionally label names
# OUTPUT: NxN correlation matrix where N is the number of ROIs
# diagnostic pdfs of motion correction and filtering results
# ----------------------------



# corr heat maps dont update w/ label names <- make sure this populates


library(methods) # this one is critical and should always be included in Rscripts

# wavelet decomposition will produce some NaNs so turn warnings off
options(warn=-1)


# get arguments
Args <- commandArgs()
self<-Args[4]
self<-substring(self,8,nchar(as.character(self)))
getPckg <- function(pckg) install.packages(pckg, repos = "http://cran.r-project.org")


# check for getopt package
pckg = try(require(getopt))

if(!pckg) 
{
	cat("Installing 'getopt' from CRAN\n")
	getPckg("getopt")
	require("getopt")
}





# ---------------------------------------------------------------------- #
# ------------------------- user input --------------------------------- #
# ......................spec function....................... #
# spec function, part of getopt library
# specifies valid input for getopt
spec = c( 
'verbose', 'v', 2, "integer" ," verbose output ",
'help'   , 'h', 0, "logical" ," print the help ", 
'rsfmri'   , 'f', 1, "character"," time-points by spatial-data csv file containing time series projections or raw data --- output of previous steps in pipeline (e.g., sub01_timeseries_mat.csv",
'dummy' , 'd', 1, "integer" ," number of dummy scans to remove",
'motion'   , 'm', 1, "character","csv file containing motion/metric nuisance variables for rsfmri data (optional)",
'nuisance'   , 'n', 1, "character","csv file containing other nuisance variables for rsfmri data (optional)",
'wavelet' , 'w', 0, "logical", "do wavelet analysis",
'bandpass'   , 'b', 1, "character"," frequences in form FLOxFHI e.g.  0.02x0.1 for bandpass filtering",
'tr'   , 't', 1, "character"," TR for the rsfMRI acquisition",
'labelnames' , 'l', 1, "character"," node label names passed in by the user, .txt format (optional)",
'output' , 'o', 1, "character"," the output prefix ")
# ...................................................... #


# .................................................... #
# reformat spec into something readable, get options
spec=matrix(spec,ncol=5,byrow=TRUE)

# get the options
opt = getopt(spec)
# .................................................... #


# ..............help message............................... #
# -h was specified
if ( !is.null(opt$help) || length(opt) == 1 )
{
	#print a friendly message and exit with a non-zero error code
	cat("\n")
	cat(paste(self,"\n"))
	for ( x in 1:nrow(spec) ) 
	{
		cat("\n")
		longopt<-paste("--",spec[x,1],sep='')
		shortopt<-paste("-",spec[x,2],sep='')
		hlist<-paste(shortopt,"|",longopt,spec[x,5],"\n \n")
		# print(hlist,quote=F)
		cat(format(hlist, width=40, justify = c("left")))
	}
	cat(format("Example: \n", width=40, justify = c("left")))
	ex<-paste(self," -f RSDATA.csv -m RSDATAMOCOparams.csv -n RSDATAcompcorr_compcorr.csv -b 0.02x0.1 -d 0 -t 2 -o OUT \n \n ")
	ex<-format(ex, width=length(ex), justify = c("left"))
	cat("\n")
	cat(ex)
	cat('Filters the timeseries data and produces an NxN correlation matrix, where N is the number of unique timeseries (the number of ROIs). The matrix is computed using either standard correlation, filtered within a specified range (the -b flag), or using wavelet analysis (-w). In the latter case, matrices are computed for each level of the analysis, where level correspondes to a particular frequency range: \n
Level 1: 0.25-0.45hz \
Level 2: 0.11-0.25hz 
Level 3: 0.06-0.11hz 
Level 4: 0.03-0.06hz 
Level 5: 0.01-0.03hz 
Level 6: 0.007-0.01hz \n\n')
	q(status=1);
}
# .................................................... #




# ....................argument check.................................# 
# check each required argument, if any are missing, quit the program

# filter range or wavelet decomposition
if ( !(is.null(opt$bandpass)) & !(is.null(opt$wavelet)) )
  {
	print(paste(" Specify EITHER wavelet decomposition or a specific frequency range to filter."))
	q(status=1)
  }



# output
if ( is.null(opt$output) ) 
  {
 	 print(paste(" Need -o option to define output. Exiting."))
  	 q(status=1)
  }



# nuisance paramters
if ( ! is.null(opt$nuisance) ) 
  {
  	print(paste("read nuis",opt$nuisance))

	# get the file extension
	ext <- sub(".*[.]", ".", opt$nuisance)
	

	# check to make sure the file is in the appropriate format
	if( (ext != ".csv") & (ext != ".txt") )
	{
		print(paste("nuisance parameters should be in .csv or .txt format"))
		print(paste("the entered file is in ", ext, " format"))
		q(status=1)
	}
  	nuisance<-read.csv(opt$nuisance)
  }



# motion parameters
if ( ! is.null(opt$motion) ) 
  {
	  print(paste("read motion",opt$motion))
	  ext <- sub(".*[.]", ".", opt$motion)
	  if( (ext != ".csv") & (ext != ".txt") )
	  {
		print(paste("motion parameters should be in .csv or .txt format"))
		print(paste("the entered file is in ", ext, " format"))
		q(status=1)
	  }

	  motion<-read.csv(opt$motion)
  }


# dummy scans
if ( ! is.null(opt$dummy) )
{
	print(paste("dummy scans ", opt$dummy))
	dummy<-opt$dummy
}




# fmri timeseries data
if ( ! is.null(opt$rsfmri) ) 
  {
	  ext <- sub(".*[.]", ".", opt$rsfmri)
	  
	  print(paste("read rsf data",opt$rsfmri))
	  if( (ext != ".csv") & (ext != ".txt") )
	  {
		print(paste("fmri data should be in .csv or .txt format, e.g. the time-series x ROI matrix."))
		print(paste("the entered file is in ", ext, " format"))
		q(status=1)
	  }

	  a<-read.csv(opt$rsfmri)
   	 
	  if ( ! is.null(opt$labelnames) ) 
   	  {
    		print(paste("rename data according to user-passed label names",opt$labelnames))
    		newnames<- read.table(opt$labelnames) 
    		colnames(a)<-newnames[,2] 
    	  }
  }
else
{
	print(paste(" Need -f option to define input rsfmri data. Exiting."))
 	 q(status=1)
}



# tr
if( is.null(opt$tr) )
{
	print(paste(" Need -t option to specify tr "))
	q(status=1)
}
# .......................................................................... #
# ------------------------------end user input-------------------------------------------- #








#.............get required packages.......................#
# ........................................................#
pckg = try(require(mFilter))

if(!pckg) 
{
	cat("Installing 'mFilter' from CRAN\n")
	getPckg("mFilter")
	require("mFilter")
}




pckg = try(require(timeSeries))

if(!pckg) 
{
	cat("Installing 'timeSeries' from CRAN\n")
	getPckg("timeSeries")
	require("timeSeries")
}




pckg = try(require(waveslim))

if(!pckg)
{
	cat("Installing 'waveslim' from CRAN\n")
	getPckg("waveslim")
	require("waveslim")
}



pckg = try(require(brainwaver))

if(!pckg)
{
	cat("Installing 'brainwaver' from CRAN\n")
	getPckg("brainwaver")
	require("brainwaver")

}
#..........................................................#
# ..................... end packages ..................... #








# ----------------------------------------------------------------------------------- #
# -------------------------- function definitions------------------------------------ #
# ----------------------------------------------------------------------------------- #


# ..................................moveFile .............................#
# simulates the bash mv function, in this case specifically moving 
# files to the diagnostic folder
# INPUT: filename
# OUTPUT: files are moved
moveFile <- function(fn)
{
	# this might be too general...
	filein <- paste("./", fn, sep="")
	fileout <- paste("./diagnostic/", fn, sep="")
	file.copy(filein, fileout)
	file.remove(filein)
}
# ........................................................................# 


# ................... filterTimeSeries .................................... #
# filter time series by a specified range, write results to pdf. used with bandpass filtering
# INPUT: time-series filtered for motion and nuisance; frequency range
# OUTPUT: filtered timeseries and diagnostic pdf
filterTimeSeries <- function(TimeSeries, freqLo, freqHi)
{
	if( level == 0)
	{
		fn<-paste(opt$output,"_filtering.pdf",sep='')
	} 
	
	# perform frequency filtering and write the data
	print(paste("Frequency filter",freqLo,"x",freqHi,"from",tr,"second TR data"))	
	voxLo<-round(1/freqLo)
	voxHi<-round(1/freqHi)
	#write.table(TimeSeries, file="TimeSeries.txt", row.names=FALSE, col.names=FALSE)
	
	# filter time series by specified range
	# only proceed if the data is valid
	fTimeSeries <- try(residuals( cffilter( TimeSeries , pl=voxHi, pu=voxLo , drift=TRUE , type="t")) )  
	if(class(fTimeSeries) == "try-error")
	{
		print( paste( "Error in filterTimeSeries"))
		return(0)
	} else
	{	
		print(paste("writing out reference filtering result",fn))
	
		pdf(fn)
		vv <- 2   # this is just picking an arbitrary ROI to display results for
		par(mfrow=c(2,2))

		# original data
		plot(TimeSeries[,vv], type='l', main='Original Time Series', ylab="Signal")
		spec.pgram( TimeSeries, taper=0, fast=FALSE, detrend=F, demean=F, log="n")
		
		# filtered data
		plot(fTimeSeries[,vv],type='l',main=paste('Filtered Time series: ',freqLo,"< f <",freqHi), ylab="Signal")
		spec.pgram(fTimeSeries, taper=0, fast=FALSE, detrend=F,demean=F, log="n")

		

		dev.off() # write pdf
		
		# move to diagnostic folder
		moveFile(fn)
		return(fTimeSeries)
	}
	
}	
# .................................................................... #


# ........................plotHeapMap....................................... #
# just a wrapper functions to make the call easier
plotHeatMap <- function( corrmat , pdfname )
{
	pdf( pdfname )
	if( is.null(opt$labelnames) )
	{
		heatmap( corrmat, symm=T, Rowv=NA )
	} else
	{
		heatmap( corrmat, symm=T, Rowv=NA, labRow=names(a), labCol=names(a) )
	}
	dev.off()

	moveFile(pdfname);
}
# ......................................................................... #



# ----------------------------------------------------------------------------------- #
# --------------------------------- end functions------------------------------------ #
# ----------------------------------------------------------------------------------- #















# ----------------------------------------------------------------------------------- #
# -------------------------------------- MAIN --------------------------------------- #
# ----------------------------------------------------------------------------------- #

#...................create plot of mocoparams.............................#
print(paste("Creating diagnostic MOCO plot..."))
moconame <- paste(opt$output, "MOCO_plot.pdf", sep='')
pdf(moconame, width=12, height=8)
par(oma=c(0,0,5,0))  # outer margins, bottom left top right

# matrix for order of plotting
m <- matrix(c(1,2,3,4), 2, 2, byrow=TRUE)
layout(m, widths=c(12,2), heights=c(5,5), TRUE)
# use:
# nf <- layout(m, widths=c(12,2), heights=c(5,5), TRUE)
# layout.show(nf)
# to see how the layout looks before plotting starts


n.tps <- dim(motion)[1]
tps <- seq(1:n.tps)   # x-axis for plotting


# 6 parameters of rigid registration
# spm/AFNI produce files with dimensions n.tps x 6
# ANTS produces a file with dimensions n.tps x 8, last 6 are rigid body motion parameters
# need to account for both options

if( dim(motion)[2] == 8)
{    motion <- motion[,3:8]   }

xtrans <- motion[,1]
ytrans <- motion[,2]
ztrans <- motion[,3]
pitch  <- motion[,4]
roll   <- motion[,5]
yaw    <- motion[,6]


# setup the range of the yaxis
ymin <- min(xtrans, ytrans, ztrans)
ymax <- max(xtrans, ytrans, ztrans)

plot(tps, xtrans, type='l', col='blue', xlab="Time Point", ylab="Translation in mm", ylim=c(ymin, ymax))
lines(tps, ytrans, type='l', col='green')
lines(tps, ztrans, type='l', col='red')

# setup the legend as a bogus plot, for proper alignment
plot(c(0,1), c(0,1), ann=F, bty='n', type='n', xaxt='n', yaxt='n')
legend("center", xpd=TRUE, legend=c("X Trans", "Y Trans", "Z Trans"), col=c('blue', 'green', 'red'), lwd=2)

# yaxis range for rotational parameters
ymin <- min(pitch, roll, yaw)
ymax <- max(pitch, roll, yaw)


# plot rotation
plot(tps, pitch, type='l', col='blue', xlab="Time Point", ylab="Rotation in mm", ylim=c(ymin, ymax))
lines(tps, roll, type='l', col='green')
lines(tps, yaw, type='l', col='red')

plot(c(0,1), c(0,1), ann=F, bty='n', type='n', xaxt='n', yaxt='n')
legend("center", xpd=TRUE, legend=c("Pitch", "Roll", "Yaw"), col=c('blue', 'green', 'red'), lwd=2)

mtext("Motion Correction Parameters", side=3, line=1, cex=2, col="black", outer=TRUE)  # figure caption

# write pdf
dev.off()

moveFile(moconame)
print(paste("Plot created"))
#........................................................................#



# ............................signal processing .............................#
# create a matrix of the time series for each ROI
# i.e., image series with 120 time points, and 4 ROIs, will produce
# a 120x4 matrix
amat<-as.matrix(a)

# so as not to mess up indexing
if ( dummy == 0 )
{
	dummy <- 1
}

inds<-c(dummy:nrow(amat)) # exclude dummy scans and reindex matrices



# even number in time series -> trigonemetric regressions only available for even n
if ( dim(amat)[1] %% 2 ==  1 ) 
	{ inds<-c((dummy+1):nrow(amat)) }


amat <- amat[inds,]
motion <- motion[inds,]
nuisance <- nuisance[inds,]




# factor out nuisance variables
# keep the portion of the signal that is not explained by nuisance vars,
# which is the residuals
if ( ! is.null(opt$nuisance) )
  {
  	print(paste("factor out nuisance vars"))
  	mp<-as.matrix(nuisance)
  	amat<-residuals(lm(amat~1+mp))
  }


# factor out motion vars
if ( ! is.null(opt$motion) )
  {
	  print(paste("factor out motion vars"))
	  mp<-as.matrix(motion)
	  mp[,1:6]<-mp[,1:6]-mp[c(1,3:nrow(mp),nrow(mp)),1:6]
	  motionconfounds<-matrix( rep(NA,2*nrow(mp)) ,nrow=nrow(mp),ncol=2)
	  motionconfounds[,1]<-sqrt( mp[,1]*mp[,1] + mp[,2]*mp[,2] + mp[,3]*mp[,3] )  # rotation magnitude
	  motionconfounds[,2]<-sqrt( mp[,4]*mp[,4] + mp[,5]*mp[,5] + mp[,6]*mp[,6] )  # translation magnitude
	  amat<-residuals(lm(amat~1+motionconfounds))
  }




# frequency filtering
print("now do frequency filtering ")
tr<-as.real(as.character(opt$tr)) 


myTimeSeries<-ts( amat ,  frequency<-1/tr )


# if the user specifies a paricular frequency range, filter for that range and save the new time series
# otherwise we will look at the different levels of the wavelet decomposition
if( ! is.null(opt$bandpass) )
{
	freqs<-c(as.numeric(unlist(strsplit(opt$bandpass,"x"))))  # splits bandpass into 2 reals
	print(paste("Frequency filter",opt$bandpass,"from",tr,"second TR data"))
	freqLo<-freqs[1]
	freqHi<-freqs[2]
	
	fTimeSeries <- filterTimeSeries(myTimeSeries, freqLo, freqHi)
	
# if bandpass filtering isn't specified, proceed to wavelet decomposition
} else 
{ 
		
	
	fTimeSeries <- myTimeSeries 
}

# write the filtered ts #
tsdf<-as.matrix(fTimeSeries)

names(tsdf)<-names(a)
#..................................................................................#



# ......................procede with analysis .............................#
# once filtering is done, create correlation matrices
if( is.null(opt$bandpass) )
{
	# use wavelet decomposition
	
	# number of time points
	n.tmPts <- dim(tsdf)[1]  


	# number of regions
	n.regions <- dim(tsdf)[2]

	# fTimeSeries is the time series matrix
	# create the wavelet correlation matrix and write to text
	print(paste("do wavelet decomposition..."))
	waveCorMat <- const.cor.list(tsdf, method = "modwt", wf = "la8", n.levels=6, boundary = "periodic", p.corr = 0.95)
	print(paste("wavelet decomposition done."))
	
	
	waveVarMat <- const.var.list(tsdf,n.levels=6)
	
	
	
	for(level in 1:6)
	{
		fn <- paste(opt$output, "_filtering_level_", level, ".pdf", sep='')
		
		name = paste(opt$output, "_cor_matrix_lvl", level, ".txt", sep="")
		
		# wavelet decomposition levels
		if(level == 1)
		{ freqLo = 0.25
		  freqHi = 0.45 } else
		if(level == 2)
		{ freqLo = 0.11
		  freqHi = 0.25 } else
		if(level == 3)
		{ freqLo = 0.06
		  freqHi = 0.11 } else
		if(level == 4)
		{ freqLo = 0.03
		  freqHi = 0.06 } else
		if(level == 5)
		{ freqLo = 0.01
		  freqHi = 0.03 } else
		if(level == 6)
		{ freqLo = 0.007 
		  freqHi = 0.01 }
	
		# read in the correlation matrix for diagnostic results
		waveCorDat <- paste("waveCorMat$d", level, sep="")
		waveCorDat <- eval(parse(text=waveCorDat))
	
		# if there are NAs in the data, its not usable
		# happens with the lower frequencies
		if( length(which(is.na(waveCorDat))) > 0)
		{
			print(paste("No valid data for level ", level))
		} else
		{
			write.table( waveCorDat , name, row.names=FALSE, col.names=FALSE)
							
			
			print(paste("writing out reference filtering result",fn))
			pdf(fn)
			vv <- 2   # this is just picking an arbitrary ROI to display results for
			par(mfrow=c(2,2))

			# plot from original, unfiltered data
			plot(tsdf[,vv], type='l', main='Original Time Series', ylab="Signal")
		
		
			spec.pgram( waveCorDat, taper=0, fast=FALSE, detrend=F, demean=F, log="n")
		
			# plot filtered data
			plot(waveCorDat[,vv],type='l',main=paste('Filtered Time series: ',freqLo,"< f <",freqHi), ylab="Signal")
			spec.pgram(waveCorDat, taper=0, fast=FALSE, detrend=F,demean=F, log="n")
			dev.off() # write pdf

			moveFile(fn)
		
			# create correlation map
			heatmapname = paste(opt$output, "_cor_matrix_map_level", level, ".pdf", sep="")
			plotHeatMap( waveCorDat, heatmapname )
		}
	}
} else
{
	# otherwise the data has alreayd been filtered so we can do a standard correlation
	im <- ( cor(tsdf) )
	
	heatmapname = paste(opt$output, "cor_matrix_map.pdf", sep="")
	name = paste(opt$output, "cor_matrix.txt", sep="")
	write.table(im, file="cor_matrix.txt", row.names=FALSE, col.names=FALSE)
	plotHeatMap( im , heatmapname)
}


#..................................................................................
# --------------------------------------------------------------------------------------- #


