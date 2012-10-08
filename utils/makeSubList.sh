#!/bin/sh

# pluta 9/19/12

# creates a text file with the path to the correlation matrix for each subject
# this file is used in group analysis scripts



usage="This script is a utility to create a text file listing the full path to each subject's correlation matrix. This is used for group analysis.

Required arguments:
	-r : the root directory for your data
	-p : the prefix for individual subjects
	-c : the name of the correlation matrix (.txt file)
	-f : the name of the output file
Optional:
	-o : subject output directory, if the correlation matrix is in a sub directory of the patient
	
	Usage: $0 -r /home/user/study/patients -p sub -o data -c wave_cor_mat_1.txt -f patlist.txt
	Output would then be:
	/home/user/study/patients/sub01/data/wave_cor_mat_1.txt
	/home/user/study/patients/sub02/data/wave_cor_mat_1.txt
	...
	
	"
	
# options
while getopts "r:p:o:c:f:h" opt ; do
  case $opt in
    r) ROOTDIR=$OPTARG ;;
    p) SUBPREFIX=$OPTARG ;;
    o) OUTPUTDIR=$OPTARG ;;
    c) CORMAT=$OPTARG ;;  
    f) OUTPUTFILE=$OPTARG ;;   
    h) echo  -e " $usage " ; exit 1 ;;
    *) echo  -e " $usage " ; exit 1 ;;
  esac
done

# make sure input is ok
if [[ ! -s $OUTPUTFILE ]]
then
	echo "No outputfile specified. Exiting."
	exit
fi

if [[ ! -s $CORMAT ]]
then
	echo "The correlation matrix has not been specified. Exiting."
	exit
fi

if [[ ! -s $ROOTDIR ]]
then
	echo "Root data directory not specified. Exiting."
fi




# want to make sure theres only one slash at the end of the path
# want explicit control of the pathname for later use in R
lastchr=${ROOTDIR#${ROOTDIR%?}}
if [[ $lastchar == "/" ]]
then
	ROOTDIR=${ROODIR%?}  # if the last character is a slash, remove it
fi

if [[ -s $OUTPUTDIR ]]
then
	lastchr=${OUTPUTDIR#${OUTPUTDIR%?}}
	if [[ $lastchar == "/" ]]
	then
		OUTPUTDIR=${OUTPUTDIR%?}
	fi
fi

# get the list of subjects
SUBLIST=`echo ${ROOTDIR}/${SUBPREFIX}*`


# create the file path here. this is the part to modify if you have an unconventional naming scheme
for SUB in $SUBLIST
do
	# for each subject, concatenate the full path
	if [[ -s $OUTPUTDIR ]]
	then
		CORMATPATH=${SUB}/${OUTPUTDIR}/${CORMAT}
	else
		CORMATPATH=${SUB}/${CORMAT}
	fi
	
	# write to output
	echo $CORMATPATH >> $OUTPUTFILE
done
