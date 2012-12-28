#!/usr/bin/env Rscript

# --------- about --------
# avants and pluta 2012
# some code is taken from the brainwaver manual, by sophie achard




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

pckg = try(require(ANTsR))
{
    cat("Requires ANTsR")
    exit 1;
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
	cat('Filters the timeseries data and produces an NxN correlation matrix, where N is the number of unique timeseries (the number of ROIs). The matrix is computed using either standard correlation, filtered within a specified range (t\n\n')
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
	  if (ext != ".nii.gz")  
	  {
		print(paste("fmri data should be a 4D .nii.gz file"))
		print(paste("the entered file is in ", ext, " format"))
		q(status=1)
	  }

	  ts<-antsImageRead(opt$rsfmri, "float", 4)
   	 
	  
  } else
{
	print(paste(" Need -f option to define input rsfmri data. Exiting."))
 	 q(status=1)
}

# brainmask
if ( ! is.null(opt$bmask) )
{
# add dimension check
# should these really be doubles?
    bmask <- antsImageRead(opt$bmask, "unsigned int", 3)
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

    fn<-paste(opt$output,"_filtering.pdf",sep='')

	
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



# ............................ main .............................#
# read in the volumes to work with:
# fmri: fmri time series (4D)
# bmask: brainmask of the average fmri image
# ROI: the seed region. create a correlation map by taking the average time-series in the seed region
# and correlating to all other voxels in the mask.

bmask <- antsImageRead( brainmask, "unsigned int", 3)
ROI <- antsImageRead( seed, "unsigned int", 3)
fmri <- antsImageRead( opt$fmri, "double", 4)

# convert fmri data to a matrix, constrained by the area under the brainmask
raw.fmri.mat <- timeseries2matrix( fmri, bmask )




# so as not to mess up indexing
if ( dummy == 0 )
{
	dummy <- 1
}

inds<-c(dummy:nrow(raw.fmri.mat)) # exclude dummy scans and reindex matrices



# even number in time series -> trigonemetric regressions only available for even n
if ( dim(raw.fmri.mat)[1] %% 2 ==  1 )
	{ inds<-c((dummy+1):nrow(raw.fmri.mat)) }

# reindex raw data excluding dummy scans
raw.fmri.mat <- raw.fmri.mat[inds,]





# factor out nuisance variables
# keep the portion of the signal that is not explained by nuisance vars,
# which is the residuals
if ( ! is.null(opt$nuisance) )
  {
    motion <- motion[inds,]
  	print(paste("factor out nuisance vars"))
  	mp<-as.matrix(nuisance)
  	fmri.mat<-residuals(lm(raw.fmri.mat~1+mp))
  }


# factor out motion vars
if ( ! is.null(opt$motion) )
  {
      nuisance <- nuisance[inds,]
	  print(paste("factor out motion vars"))
	  mp<-as.matrix(motion)
	  mp[,1:6]<-mp[,1:6]-mp[c(1,3:nrow(mp),nrow(mp)),1:6]
	  motionconfounds<-matrix( rep(NA,2*nrow(mp)) ,nrow=nrow(mp),ncol=2)
	  motionconfounds[,1]<-sqrt( mp[,1]*mp[,1] + mp[,2]*mp[,2] + mp[,3]*mp[,3] )  # rotation magnitude
	  motionconfounds[,2]<-sqrt( mp[,4]*mp[,4] + mp[,5]*mp[,5] + mp[,6]*mp[,6] )  # translation magnitude
	  fmri.mat<-residuals(lm(fmri.mat~1+motionconfounds))
  }




# frequency filtering
print("now do frequency filtering ")
tr<-as.real(as.character(opt$tr)) 


myTimeSeries<-ts( fmri.mat ,  frequency<-1/tr )


freqs<-c(as.numeric(unlist(strsplit(opt$bandpass,"x"))))  # splits bandpass into 2 reals
	
freqLo<-freqs[1]
freqHi<-freqs[2]
	
fTimeSeries <- filterTimeSeries(myTimeSeries, freqLo, freqHi)
	

main.ts <- as.matrix(fTimeSeries)

bmask.range <- (bmask == 1)
roi.range <- (ROI == 1)

roi.ts <- fTimeSeries[roi.range]
roi.ts <- rowMeans(roi.ts)

ntwk <- rep( NA, ncol(main.ts) )

for( i in 1:ncol(main.ts) )
{
    corr <- cor.test(roi.ts, main.ts[,i])
    ntwk[i] <- corr$est
}

corr.img <- array(0, dim=c(xdim, ydim, zdim) )
corr.img[bmask.range] <- as.real(ntwk)
corr.img <- as.antsImage(corr.img)
antsImageWrite(corr.img, name)




write.table(tsdf, file="filtered_timeseries.txt")
# but, for wavelet decomposition, this wont be truly filtered...
# this is the filtered time series of the ROI, remember
#..................................................................................#



#..................................................................................
# --------------------------------------------------------------------------------------- #


