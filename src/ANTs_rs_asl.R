
/mnt/data/avants/Wolk/fmri/ipmi2012

using restingstate_network.R & newAntsExample.sh

########################################################
# network analysis with asl and bold-from-asl = B.Avants Oct 2012
########################################################
# draw the roi in label 2 on the mask
mask<-antsImageRead( paste(pre,"cortexb.nii.gz",sep=''), 'float', 3 )
wb<-( mask > 0 )           # whole brain
nwb<-( mask < 1 )         # not brain
wbvec<-mask[ wb ]       # mask data in whole brain
roivec<-( wbvec > 1  )    # where the ROI exists
wbvec[]<-1                    # redefine whole brain as = 1 , vs >= 1
saslimg<-antsImageRead( fn, 'float', 4 )
SmoothImage( 4,  asl, 1, saslimg) # just a little smoothing in space/time
gmat <- timeseries2matrix( asl , wb  )  # make a matrix of asl (not smoothed)
#  below , we assume we already ran motion-corr
motionparamsandcompcorr<-predictors$nuis[,3:8] # motion+compcorr
gmatasl <-residuals( lm( gmat ~ motionparamsandcompcorr ))
gmatbold<-residuals( lm( gmat ~ predictors$xideal + motionparamsandcompcorr ) )
leftinds<-shift(c(1:nrow(gmat)),1)
rightinds<-shift(c(1:nrow(gmat)),-1)
# surround subtraction - difference of me with my neighbors' average
gmatasl<-gmatasl - 0.5 * ( gmatasl[leftinds,] + gmatasl[rightinds,] )
taginds<-c(1:(nrow(gmat)/2))*2
controlinds<-taginds-1
gmatasl[controlinds,]<-gmatasl[controlinds,]*(-1)
# ok! done w/asl specific stuff , now do standard network analysis

let me know if you have questions.

brian
