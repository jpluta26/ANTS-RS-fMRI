# Processing steps
# 
# see help in individual scripts for more detail
#
#
# 0) Preprocess data (using SPM, pipedream, etc)
# generate binary brainmask, 3 tissue segmentations, and ROIs in T1 space
# generate a 4D-nifti file of functional data
#
# 1) ANTs_fmriToT1.sh
# motion correction of fmri data, coregistration between T1 and fmri
# put everythign into fmri space
#
#
# $>~/scripts/ANTs_fmriToT1.sh -i fmri_series.nii.gz -o sub01 -t sub01_T1.nii.gz -r sub01_ROI.nii.gz -c sub01_csf.nii.gz -w sub01_wm.nii.gz -g sub01_gm.nii.gz -b sub01_brainmask.nii.gz

# 
#
# 2) ANTs_rsfmri_pre.sh
# extract the time series under each ROI
#
# $>~/scripts/ANTs_rsfmri_pre.sh -i fmri_series_MOCO.nii.gz -o sub01 -r sub01_ROI_fmri.nii.gz -c sub01_csf_fmri.nii.gz -g sub01_gm_fmri.nii.gz -w sub01_wm_fmri.nii.gz -b sub01_brainmask_fmri.nii.gz -a sub01avg.nii.gz
#
#
#
# 3) ANTs_rsfmri_filter.sh
# filter each time series for motion, nuisance variables (csf, wm, global signal), and filter by a specified frequency. create a correlation matrix (or matrices) using correlation or wavelet decomposition. generate SPM-style plot of motion correction parameters
#
# $>~/scripts/ANTs_rsfmri_filter.sh -f sub01.csv -d 0 -m sub01MOCOparams.csv -n sub01compcorr_compcorr.csv -tr 3 -w -o sub01

