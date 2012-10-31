# this script uses antsr to combine and antsrsfmripre and antsrsfmrifilter into one script

#!/usr/bin/env Rscript



# check for required packages
pckg = try(require(getopt))

if(!pckg)
{
    cat("Installing 'getopt' from CRAN\n")
    install.packages("getopt", repos="http://cran.r-project.org")
}

pckg = try(require(ANTsR))

if(!pckg)
{
    # installing ANTsR isn't as simple, the user needs to do this
    cat("ANTsR is required. See www.instructions.com for install instructions.")
    exit 1;
}



# load required libraries
library(methods)
library(getopt)
library(ANTsR)

# --------------------- user input ------------------------ #

# .................... spec function ...................... #
# specify input for getopt
spec = c(
'fmri', 'f', 1, "character", " motion-corrected functional image series (4D) ",
'roi', 'r', 1, "character", " an integer-valued 3D file containing the ROIs ",
'csf', 'c', 1, "character", " csf probability map ",
'graymatter', 'c', 1, "character", " graymatter probability map ",
'whitematter', 'c', 1, "character", " whitematter probability map",
'output', 'o', 1, "character", " output prefix",
'wavelet', 'w', 0, "logical", " do wavelet analysis",
'dummy', 'd', 1, "integer", " number of dummy scans to exclude (default=0)",
'bandpass', 'b', 1, "character", " frequency range to filter data, in the form LOxHI, e.g. 0.02x0.1 for bandpass filtering",
'labelnames', 'l', 1, "character", " ROI label names in .txt format (optional)")

# nuisance, motion will be created
# maybe TR can be taken from the image data
# ........................................................... #


# .................................................... #
# reformat spec into something readable, get options
spec=matrix(spec,ncol=5,byrow=TRUE)

# get the options
opt = getopt(spec)
# .................................................... #

### update this!
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


# ................ argument check .................... #
# make sure required arguments are supplied and valid

