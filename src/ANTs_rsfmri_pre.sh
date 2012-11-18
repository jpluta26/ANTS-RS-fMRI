#!/bin/bash

# by brian avants and john pluta 6/11/12

# script to extract the average time series for each ROI. that is, take the average value of each each time point, under
# each ROI. the time series are stored in a TxN matrix text file, where T is the number of time points and N is the number of
# ROIs. optionally performs motion-correction, if this hasn't been performed in a prior step.

# INPUT: brainmask, 3 tissue probability maps, ROI, all in fmri space. optionally, average fmri image.
# OUTPUT: TxN text file of time-series data. diagnostic images if motion correction is performed.





# necessary paths - user defined

# at CFN use:
#export ANTSPATH=/usr/local/ANTS/ANTs-1.9.x-Linux/bin/:$ANTSPATH
#export LD_LIBRARY_PATH=/usr/local/ANTS/ANTs-1.9.x-Linux/bin:$LD_LIBRARY_PATH

# at PICSL use:
export ANTSPATH=/home/avants/bin/ants/:${ANTSPATH}



# ------------------------------ usage and input -------------------------------- #
# check and make sure all the necessary programs are available
for prog in antsMotionCorr  ImageMath N3BiasFieldCorrection Atropos ThresholdImage sccan ; do
  command -v $prog >/dev/null 2>&1 || { echo -e >&2 "This program $0 requires $prog but it's not installed. \n please add the ANTS binary directory to your path.  \n see brianavants.wordpress.com to find ANTs installation details .  Aborting."; exit 1; }
done
usage=" \n $0 takes an fmri series as input, motion corrects it, 
        segments it and outputs it to a csv file named  outputdir/outputprefix.csv. 
        alternatively , if you pass in an ROI image, the ROI image will be masked 
        with the rsfmri-derived gray matter image and the output will be the average 
        time series in each roi --- written to outputdir/outputprefix.csv  \n 

	required input:
	-i : the fmri time series (4D image)
	-o : output prefix
	-r : the ROIs for which connectivity is computed. all ROIs should be contained in a single file, and each should have a unique integer value (3D image)	
	-c : CSF prior        
	-g : GM prior 	
	-w : WM prior 
	-b : brainmask, should be a binary image
	
	optional input:
	-a : the average of the fmri time series (3D). use this option if you have already performed motion correction.
	
	All images should be in the same T2 space.
	
	
	All image files should be in .nii.gz format

	outputs have extensions
        
       MOCOparams.csv --- motion correction and metric measurements 
     corrected.nii.gz --- motion corrected image (4d)
           avg.nii.gz --- rsfmri average image (3d)
         bmask.nii.gz --- brain mask
           seg.nii.gz --- 3 tissue segmentation
         *compcorr*   --- compcorr related outputs including compcorred image
         *rsfmri*jpg  --- diagnostic images 

        Usage: $0 -i inputimage.nii.gz -r roi.nii.gz -b brainmask -c csfpriorimage -g gmpriorimage -w wmpriorimage -o outputdir/outputprefix \n \n
        inputimage is a 4D image 
        roi.nii.gz is a 3D integer valued label image \n \n "

# .................. get user input ................... #
while getopts "i:o:r:c:g:w:b:t:a:h" opt ; do
  case $opt in
    i) img=$OPTARG ;;
    o) out=$OPTARG ;;
    r) roi=$OPTARG ;;
    c) CSFPRIOR=$OPTARG ;;  # prob01
    g) GMPRIOR=$OPTARG ;;   # prob02
    w) WMPRIOR=$OPTARG ;;   # prob03
    b) BMASK=$OPTARG ;; 
	a) AVG=$OPTARG ;;
    h) echo  -e " $usage " ; exit 1 ;;
    *) echo  -e " $usage " ; exit 1 ;;
  esac
done
# ..................................................... #

# ....................checkDim ........................ #
function checkDim
{
    #  the sole input is the image to check dimensions; check against the dimensions of the average fmri image	
    IMG=$1

    str=`PrintHeader $IMG | grep 'Voxel Spacing' | awk '{print $4, $5, $6}' | tr '[' ' ' | tr ']' ' '`
    IMGXDIM=`echo $str | awk '{print $1}' | tr ',' ' '`
    IMGYDIM=`echo $str | awk '{print $2}' | tr ',' ' '`
    IMGZDIM=`echo $str | awk '{print $3}'`

    if [ $IMGXDIM != $T2XDIM ] 
    then
	echo " "
        echo "$img and $IMG must have the same dimensions - see usage. Exiting"
	echo " "
    	exit 1
    elif [ $IMGYDIM != $T2YDIM ] 
    then
   	echo " "
	echo "$img and $IMG must have the same dimensions - see usage. Exiting"
    	echo " "
	exit 1
    elif [ $IMGZDIM != $T2ZDIM ] 
    then
        echo " "
	echo "$img and $IMG must have the same dimensions - see usage. Exiting"
    	echo " "
	exit 1
    fi
}
# ..................................................... #

# ................ check for required arguments ............. #
if [[ ! -s $img ]] 
then
  	echo " "
	echo "-e the input image $img does not exist, exiting ... pass -h option on command line"
	echo " "
	exit 1
else
    str=`PrintHeader $img | grep 'Voxel Spacing' | awk '{print $4, $5, $6}' | tr '[' ' ' | tr ']' ' '`
    T2XDIM=`echo $str | awk '{print $1}' | tr ',' ' '`
    T2YDIM=`echo $str | awk '{print $2}' | tr ',' ' '`
    T2ZDIM=`echo $str | awk '{print $3}' | tr ',' ' '`
fi

if [[ ! -s $CSFPRIOR ]]
then
	echo " "
	echo "-e the CSF prior $CSFPRIOR does not exist, exiting ... pass -h option on command line"
	echo " "
	exit 1
else
    checkDim $CSFPRIOR
fi

if [[ ! -s $GMPRIOR ]] 
then
	echo " "
	echo "-e the GM prior $GMPRIOR does not exist, exiting ... pass -h option on command line"
	echo " "
	exit 1
else
    checkDim $GMPRIOR
fi

if [[ ! -s $WMPRIOR ]] 
then
	echo " "
	echo "-e the WM prior $WMPRIOR does not exist, exiting ... pass -h option on command line"
	echo " "
	exit 1
else
    checkDim $WMPRIOR
fi

if [[ ! -s $BMASK ]]
then
	echo " "
	echo "-e the brainmask $BMASK does not exist, exiting ... pass -h option on command line"
	echo " "
	exit 1
else
    checkDim $BMASK
fi

# for now, ROI is required
if [[ ! -s $roi ]]
then
  	echo " "
  	echo "the input image $roi does not exist, exiting ... pass -h option on command line"
	echo " "
  	exit 1
else
    checkDim $roi
fi



if [[ ${#out} -lt 1 ]] 
then
	echo " "
	echo "no output defined  ... pass -h option on command line"
  	echo -e  $usage
	echo " "
  	exit 1
fi

 
bd=$(dirname $out) 


if [[ ! -s $bd ]]
then
	echo " "
	echo "the output directory $bd does not exist, creating it"
	echo " "
  	mkdir -p $bd 
fi
# .................................................................... #



# ----------------------------- main ---------------------------------- #
echo img is $img basedir is $bd outprefix is $out 

# if an average image isn't specified, do motion correction and produce diagnostics
if [[ ! -s $AVG ]]
then
  echo " "
  echo "Performing motion correction on $img "
  echo " "
  
  antsMotionCorr  -d 3 -a $img -o ${out}avg.nii.gz  #fortex
  ExtractSliceFromImage 4 $img ${out}_slice_original.nii.gz 0 32
  antsMotionCorr  -d 3 -o [${out},${out}_MOCO.nii.gz,${out}avg.nii.gz] -m mi[ ${out}avg.nii.gz , $img , 1 , 32 , Regular, 0.01 ]  -t Rigid[ 0.25 ] -i 50 -u 1 -e 1 -s 1 -f 1 -n 10 -l 1 #fortex
  ExtractSliceFromImage 4 ${out}_MOCO.nii.gz ${out}_slice_corrected.nii.gz 0 32
  ExtractSliceFromImage 3 ${out}_slice_corrected.nii.gz  ${out}_slice_corrected_x.nii.gz 1 22
  ExtractSliceFromImage 3 ${out}_slice_original.nii.gz  ${out}_slice_original_x.nii.gz 1 22
  ConvertToJpg ${out}_slice_original_x.nii.gz ${out}rsfmri_discon.jpg 
  ConvertToJpg ${out}_slice_corrected_x.nii.gz ${out}rsfmri_smooth.jpg 


  DIAG=${outdir}/diagnostic
  mkdir $DIAG
  mv *.jpg $DIAG
  mv *corrected_x.nii.gz $DIAG
  mv *original_x.nii.gz $DIAG
  mv *slice_corrected.nii.gz $DIAG
  mv *slice_original.nii.gz $DIAG
else cp $AVG ${out}avg.nii.gz
fi
	
ThresholdImage 3 $BMASK ${out}bmask.nii.gz 0.5 1.e9
  
# rename tissue priors for ease of use  
cp $CSFPRIOR ${out}prob01.nii.gz
cp $GMPRIOR ${out}prob02.nii.gz
cp $WMPRIOR ${out}prob03.nii.gz



ImageMath 3 ${out}bmask.nii.gz GetLargestComponent ${out}bmask.nii.gz #fortex
ImageMath 4 ${out}compcorr.nii.gz CompCorrAuto ${out}.nii.gz ${out}bmask.nii.gz 2 #fortex

# bias correction
if [[ ! -s  ${out}n3.nii.gz ]] 
then   
	N3BiasFieldCorrection 3 ${out}avg.nii.gz ${out}n3.nii.gz 2 ${out}bmask.nii.gz ; 
fi

# 3 tissue segmentation of fmri image, guided by input tissue priors and brainmask
echo " "
echo "performing tissue segmentation of fMRI image..."
echo " "
Atropos -d 3 -a ${out}n3.nii.gz -a ${out}compcorr_variance.nii.gz -m [0.1,1x1x1] -o [${out}seg.nii.gz,${out}prob%02d.nii.gz] -x ${out}bmask.nii.gz -c [10,0] -i priorprobabilityimages[3,${out}prob%02d.nii.gz,0.5]
echo "segmentation ok"

# if you have an ROI file, extract time series from each label
echo " "
echo "extracting time series from $roi"
echo " "
 
if [[ -s $roi ]] ; then
  ThresholdImage 3 ${out}prob02.nii.gz ${out}temp.nii.gz 0.1 1   # restrict analysis to cortical regions
  MultiplyImages 3 ${out}temp.nii.gz $roi ${out}roi.nii.gz 
  rm  ${out}temp.nii.gz
  sccan --timeseriesimage-to-matrix [ ${out}_MOCO.nii.gz , ${out}roi.nii.gz , 0 , 0 ] -o ${out}_timeseries_mat.csv 
  exit 
fi

# if no ROI is specified, use a data driven approach
# at this point an ROI is required- this approach is too memory intensive and will cause a crash
# if [[ ! -s ${out}.csv ]] ; then 
#  sccan --timeseriesimage-to-matrix [ ${out}.nii.gz , ${out}bmask.nii.gz , 0 , 0 ] -o ${out}.csv #fortex
#  exit 
#  sccan --svd cgspca[ ${out}.csv ,  ${out}bmask.nii.gz , -0.05 ] -n 100 -i 1 --PClusterThresh 20 -o ${out}RSFNodes.nii.gz  #fortex
# fi
# exit 

# generates the label set generated by spca 
  MultiplyImages 3 ${out}avg.nii.gz 0 ${out}temp.nii.gz
  ct=1
  for x in ${out}RSFNodesView1vec*gz ; do
    ThresholdImage 3 $x ${out}temp2.nii.gz 1.e-6 999
    MultiplyImages 3 ${out}temp2.nii.gz $ct ${out}temp2.nii.gz
    ImageMath 3 ${out}temp.nii.gz addtozero ${out}temp.nii.gz ${out}temp2.nii.gz
    let ct=$ct+1
  done
  cp ${out}temp.nii.gz ${out}nodes.nii.gz
  rm ${out}temp.nii.gz ${out}temp2.nii.gz
  rm ${out}n3.nii.gz
  rm ${out}bmask.nii.gz
  rm ${out}roi.nii.gz
exit
# ------------------------------------------------------------------------- #
