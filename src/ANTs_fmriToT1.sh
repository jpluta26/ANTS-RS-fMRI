#!/bin/bash

# john pluta and brian avants 6/19/2012

# script to handle T1/T2 coregistration. the script assumes the T1 mprage image has been processed using some other method 
# (pipedream, spm, fsl, etc).
# INPUT: binary brainmask, 3 tissue probability maps, ROI file, all in T1 space. functional image series as one 4-dimensional file.
# OUTPUT: brainmask, probability maps, and ROI in functional image space. motion-corrected functional image series (4D) and an average image
# of all the functional images (3D). images will be appended with "_T2".


# TODO: remove dependency on c3d - can probably replace this with PrintHeader
# remove bet2 dependency



# necessary paths - user defined
export ANTSPATH=/usr/local/ANTS/ANTs-1.9.x-Linux/bin/:$ANTSPATH
export LD_LIBRARY_PATH=/usr/local/ANTS/ANTs-1.9.x-Linux/bin:$LD_LIBRARY_PATH




# ----------------------- usage and input ---------------------------- #
# check and make sure all the necessary programs are available
for prog in antsMotionCorr bet2 ImageMath N3BiasFieldCorrection Atropos ThresholdImage sccan ; do
  command -v $prog >/dev/null 2>&1 || { echo -e >&2 "This program $0 requires $prog but it's not installed. \n please add the ANTS binary directory to your path.  \n see brianavants.wordpress.com to find ANTs installation details .  Aborting."; exit 1; }
done
usage=" \n This script creates the necessary files for FC preprocessing from the space of the T1 image.

	   -Motion corrects time series and generates average functional image
	   -Coregisters average functional image to T1
	   -Transforms brainmask, tissue priors, and ROIs to functional image space  \n 

        required input:
	-i : the fmri time series (4D image)
	-o : output prefix

	the following arguments are required and should be in T1 space:
	-t : the T1 image	
	-r : the ROIs for which connectivity is computed. all ROIs should be contained in a single file, and each should have a unique integer value (3D image)	
	-c : CSF prior        
	-w : WM prior
	-g : GM prior
	-b : brainmask, should be a binary image

	
	
	All image files should be in .nii.gz format

	Usage: $0 -i 4DtimeSeries.nii.gz -t T1img.nii.gz -r T1roi.nii.gz -b T1brainmask -c T1csfpriorimage -g T1gmpriorimage -w T1wmpriorimage -o outputdir/outputprefix \n \n"
        
# ............... get user input ............. #
while getopts "i:o:r:c:g:w:b:t:h" opt ; do
  case $opt in
    i) img=$OPTARG ;;
    o) out=$OPTARG ;;
    r) roi=$OPTARG ;;
    c) CSFPRIOR=$OPTARG ;;  # prob01
    g) GMPRIOR=$OPTARG ;;   # prob02
    w) WMPRIOR=$OPTARG ;;   # prob03
    b) BMASK=$OPTARG ;; 
    t) T1=$OPTARG ;;
    h) echo  -e " $usage " ; exit 1 ;;
    *) echo  -e " $usage " ; exit 1 ;;
  esac
done
# ............................................. #

# ........  check that the necessary arguments exist ...... #
if [[ ! -s $img ]] ; then
  echo -e the input image $img does not exist, exiting ... pass -h option on command line
  exit 1
fi


if [[ ! -s $roi ]] ; then
  echo the input image $roi does not exist
fi


if [[ ${#out} -lt 1 ]] ; then
  echo no output defined  ... pass -h option on command line
  echo -e  $usage
  exit 1
fi

 
if [[ ! -s $CSFPRIOR ]] ; then
	echo CSF image does not exist, exiting ... pass -h option on command line
	exit 1
fi


if [[ ! -s $GMPRIOR ]] ; then
	echo GM image does not exist, exiting ... pass -h option on command line
	exit 1
fi


if [[ ! -s $WMPRIOR ]] ; then
	echo WM image does not exist, exiting ... pass -h option on command line
	exit 1
fi

if [[ ! -s $BMASK ]] ; then
	echo brainmask does not exist, exiting ... pass -h option on command line
	exit 1
fi

if [[ ! -s $T1 ]] ; then
	echo T1 image does not exist, exiting ... pass -h option on command line
	exit 1
fi


outdir=$(dirname $out) 


if [[ ! -s $outdir ]] ; then
  echo the output directory $outdir does not exist, creating it
  mkdir -p $outdir 
fi
# ............................................................

# ------------------------------------------------------------------- #



# ----------------------------- main ------------------------------ #
echo img is $img basedir is $outdir outprefix is $out 

AVG=${out}avg.nii.gz
echo "AVG is $AVG"


# motion correction and diagnostics
if [[ ! -s $AVG ]] ; then
  echo "Performing motion correction on $img "
  antsMotionCorr  -d 3 -a $img -o ${out}avg.nii.gz  #fortex
  ExtractSliceFromImage 4 $img ${out}_slice_original.nii.gz 0 32
  antsMotionCorr  -d 3 -o [${out},${out}.nii.gz,${out}avg.nii.gz] -m mi[ ${out}avg.nii.gz , $img , 1 , 32 , Regular, 0.01 ]  -t Rigid[ 0.25 ] -i 50 -u 1 -e 1 -s 1 -f 1 -n 10 -l 1 #fortex
  ExtractSliceFromImage 4 ${out}.nii.gz ${out}_slice_corrected.nii.gz 0 32
  ExtractSliceFromImage 3 ${out}_slice_corrected.nii.gz  ${out}_slice_corrected_x.nii.gz 1 22
  ExtractSliceFromImage 3 ${out}_slice_original.nii.gz  ${out}_slice_original_x.nii.gz 1 22
  ConvertToJpg ${out}_slice_original_x.nii.gz ${out}rsfmri_discon.jpg 
  ConvertToJpg ${out}_slice_corrected_x.nii.gz ${out}rsfmri_smooth.jpg 

fi



# ................ coregistration ..........................
# get image dimensions
str=`c3d ${out}avg.nii.gz -info-full | grep 'Voxel Spacing' | awk '{print $4, $5, $6}' | tr '[' ' ' | tr ']' ' '`
XDIM=`echo $str | awk '{print $1}' | tr ',' ' '`
YDIM=`echo $str | awk '{print $2}' | tr ',' ' '`
ZDIM=`echo $str | awk '{print $3}'`

# skull strip T1 and T2 for coregistration
# coregistration is better if we skull strip the average fmri image
ImageMath 3 ${out}tempT1.nii.gz m $T1 $BMASK
bet2 ${out}avg ${out}tempT2 -t

# create a lo-res image for coregistration
ResampleImageBySpacing 3 ${out}tempT1.nii.gz ${out}t1_lores.nii.gz $XDIM $YDIM $ZDIM


FIXED=${out}t1_lores.nii.gz
MOVING=${out}tempT2.nii.gz



# perform coregistration 
exe="antsRegistration -d 3 -r [ $FIXED , $MOVING, 1] -o [ ${out}toT1 , ${out}toT1.nii.gz ] -m Mattes[ $FIXED, $MOVING, 1, 32, Regular, 0.5] -t Translation[0.5] -f 6x4x2 -s 3x2x1 -c [100x100x100, 1e-08, 10]"
echo $exe
$exe

# transform files to fmri space
count=0
for IMG in $BMASK $CSFPRIOR $GMPRIOR $WMPRIOR $roi
do
	# setup the image names
	# if the images are in another directory, we need to parse the text
	# just keep the image name and prepend ${out}
	arr=(`echo $IMG | tr "/" " "`)	  # parse text into an array
	len=${#arr[*]}	          # get the length of the array
	len=$((len-1))	                  # adjust the index
	NAME=`basename ${arr[$len]} .nii.gz`
	
	# make sure we don't duplicate sub id if that is also the output name
	# ie don't want sub01_brainmask.nii.gz to be written as sub01sub01_brainmask_T2.nii.gz

	# if $out exists in $IMG, we don't prepend it.
	if echo "$IMG" | grep -q "$out"; 
	then	
		OUT=${NAME}_T2.nii.gz 
	else	
		OUT=${out}${NAME}_T2.nii.gz
	fi
	

	count=$((count+1))
	
	# apply transformation
	cmd="antsApplyTransforms -d 3 -i $IMG -r ${out}avg.nii.gz -o $OUT -t [${out}toT11Translation.mat , 1] "

	# for ROI file
	if [[ $count -eq 5 ]] ;
	then
		cmd=${cmd}" -n NearestNeighbor"
	fi

	$cmd
done
# ............................................................. #


# .................... quality control ........................ #
#
# do a quick quality assessment on the registration; 
# binarize the average image and make sure there is overlap between the binary
# image and the brainmask. this will not give particularly specific information
# on the quality of the registration, but it can detect a failure.
BMASK=${out}_brainmask_T2.nii.gz

# this may change depending on the intensity range of your images
# 100 is good for BOLD

# ? maybe i can extract this from the image histogram
THRESHLO=100

NAME=`basename $AVG .nii.gz`
AVGIMGBIN=${NAME}bin.nii.gz

c3d $AVG -thresh $THRESHLO inf 1 0 -o $AVGIMGBIN

# get dice overlap, assign to OVL
OVL=$(c3d $AVGIMGBIN $BMASK -overlap 1 | awk '{print $6}' )

# remove the comma
OVL=`echo "${OVL%?}"`

OVLTHRESH=0.8  # dice overlap threshold. if overlap between the binary average image and the brainmask 
		# is less than the threshold, theres been an error in the registration process.

# if OVL > OVLTHRESH, count=0
# if OVL < OVLTHRESH, count=1
dif=`echo $OVL - $OVLTHRESH | bc`
count=`echo $dif | grep "-" | wc -l`

if [[ $count -eq 0 ]]
then
	echo "Coregistration ok."
else
	echo "Warning!!!! Coregistration failed."
fi



# clean up
rm ${out}t1_lores.nii.gz
rm ${out}tempT1.nii.gz
rm ${out}tempT2.nii.gz
rm $AVGIMGBIN
# .................................................... #

# --------------------------------------------------------------- #