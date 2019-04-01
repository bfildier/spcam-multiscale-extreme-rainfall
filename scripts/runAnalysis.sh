#!/bin/bash

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $scriptdir

# Load modules if necessary
source $scriptdir/loadPlatformSpecific.sh analysis

#-----------------------------#
#   Compute all statistics    #
#-----------------------------#

if [[ "$HOSTNAME" == "jollyjumper" ]]
then
	dataroot="/Users/bfildier/Data"
elif [[ "$HOSTNAME" == "edison"* ]] || [[ "$HOSTNAME" == "cori"* ]]
then	
	dataroot=${scriptdir%\/*}
fi

# compsets="FSPCAMm_AMIP FAMIPC5"
compsets="FSPCAMm_AMIP"
experiments="piControl abrupt4xCO2"
# subsets="tropics ocean land"
subsets="tropics"
pr_ids="PRECT CRM_PREC_I90 CRM_PREC_I75 CRM_PREC_I50 CRM_PREC_I25 CRM_PREC_I10"
pr_ids_2D="PRECT CRM_PREC_I90 CRM_PREC_I75 CRM_PREC_I50"
percentiles="90.0      92.0567  93.6904  94.9881  96.0189  96.8377  97.4881  98.0047"\
"  98.4151  98.7411  99.0      99.2057  99.369   99.4988  99.6019  99.6838"\
"  99.7488  99.8005  99.8415  99.8741  99.9     99.9206  99.9369  99.9499"\
"  99.9602  99.9684  99.9749  99.98    99.9842  99.9874  99.99    99.9921"\
"  99.9937  99.995   99.996   99.9968  99.9975  99.998   99.9984  99.9987"
var2D_ids="CAPE CBMF CIN CRM_PREC_I10 CRM_PREC_I25 CRM_PREC_I50 CRM_PREC_I75 CRM_PREC_I90"\
" ICEFRAC MF PCLDBOT PCLDTOP PRECC PRECL PRECSH PRECT PRECTMX PS SST TMQ TS cin_Cu cinlcl_Cu"\
" PRECFRAC_I90 PRECFRAC_I75 PRECFRAC_I50 PRECFRAC_I25 PRECFRAC_I10"\
" PRECAREA_I90 PRECAREA_I75 PRECAREA_I50 PRECAREA_I25 PRECAREA_I10"
# var3D_ids="CRM_QC_I10 CRM_QC_I25 CRM_QC_I50 CRM_QC_I75 CRM_QC_I90"\
# " CRM_QPC_I10 CRM_QPC_I25 CRM_QPC_I50 CRM_QPC_I75 CRM_QPC_I90"\
# " CRM_T_I10 CRM_T_I25 CRM_T_I50 CRM_T_I75 CRM_T_I90"\
# " CRM_W_I10 CRM_W_I25 CRM_W_I50 CRM_W_I75 CRM_W_I90"\
# " EVAPPREC FU FV OMEGA PDELDRY Q RELHUM SPMC SPMCUP T U V"
var3D_ids="OMEGA_CRM_W_I90 RHO_CRM_T_I90 CRM_OMEGA_WT_I90 QVSATENV_CRM_T_I90"\
" OMEGA_CRM_W_I75 RHO_CRM_T_I75 CRM_OMEGA_WT_I75 QVSATENV_CRM_T_I75"\
" OMEGA_CRM_W_I25 RHO_CRM_T_I25 CRM_OMEGA_WT_I25 QVSATENV_CRM_T_I25"\
" OMEGA_CRM_W_I10 RHO_CRM_T_I10 CRM_OMEGA_WT_I10 QVSATENV_CRM_T_I10"


## Compute precipitation statistics for all prepcipitation variables
computeAllPrecipQuantiles=false
if [ "$computeAllPrecipQuantiles" == "true" ]
then
	for compset in `echo $compsets`
	do
		for experiment in `echo $experiments`
		do
			for subset in `echo $subsets`
			do
				echo 
				echo "------------------------------------------------------------------------"
				echo "     Computing precipitation statistics for:"
				echo "     - compset $compset"
				echo "     - experiment $experiment"
				echo "     - subset $subset"
				echo "     - precipitation variables ${pr_ids}"
				echo "------------------------------------------------------------------------"
				echo
				python computePrecipitationQuantiles.py ${pr_ids} -d ${dataroot} -c $compset -e $experiment -s $subset
			done
		done
	done
fi

## Compute joint precipitation statistics for all prepcipitation variables
computeJointPrecipStats=false
if [ "$computeJointPrecipStats" == "true" ]
then
	for compset in `echo $compsets`
	do
		for experiment in `echo $experiments`
		do
			for subset in `echo $subsets`
			do
				echo 
				echo "------------------------------------------------------------------------"
				echo "     Computing joint precipitation statistics for:"
				echo "     - compset $compset"
				echo "     - experiment $experiment"
				echo "     - subset $subset"
				echo "     - precipitation variables ${pr_ids}"
				echo "------------------------------------------------------------------------"
				echo
				python computeJointPrecipitationStatistics.py -p1 ${pr_ids_2D} -p2 ${pr_ids_2D} \
				 -d ${dataroot} -c $compset -e $experiment -s $subset
			done
		done
	done
fi

## Compute means for all 2D variables
computeAllMeans=false
if [ "$computeAllMeans" == "true" ]
then
	for compset in `echo $compsets`
	do
		for experiment in `echo $experiments`
		do
			for subset in `echo $subsets`
			do
				echo 
				echo "------------------------------------------------------------------------"
				echo "     Computing means of 2D variables for:"
				echo "     - compset $compset"
				echo "     - experiment $experiment"
				echo "     - subset $subset"
				echo "     - target variables ${var2D_ids}"
				echo "------------------------------------------------------------------------"
				echo
				python computeMean.py ${var2D_ids} -d ${dataroot} -c $compset -e $experiment -s $subset
			done
		done
	done
fi

## Compute mean profiles for all 3D variables
computeAllMeanProfiles=false
if [ "$computeAllMeanProfiles" == "true" ]
then
	for compset in `echo $compsets`
	do
		for experiment in `echo $experiments`
		do
			for subset in `echo $subsets`
			do
				echo 
				echo "------------------------------------------------------------------------"
				echo "     Computing mean profiles of 3D variables for:"
				echo "     - compset $compset"
				echo "     - experiment $experiment"
				echo "     - subset $subset"
				echo "     - target variables ${var3D_ids}"
				echo "------------------------------------------------------------------------"
				echo
				python computeMeanProfile.py ${var3D_ids} -d ${dataroot} -c $compset -e $experiment \
				 -s $subset
			done
		done
	done
fi

## Compute mean at locations of precipitation quantile
computeAllMeansAtQuantile=false
if [ "$computeAllMeansAtQuantile" == "true" ]
then
	for compset in `echo $compsets`
	do
		for experiment in `echo $experiments`
		do
			for subset in `echo $subsets`
			do
				echo 
				echo "------------------------------------------------------------------------"
				echo "     Computing means of 2D variables for:"
				echo "     - compset $compset"
				echo "     - experiment $experiment"
				echo "     - subset $subset"
				echo "     - target variables ${var2D_ids}"
				echo "     over location of percentiles $percentiles"
				echo "     for precipitation variables ${pr_ids}"
				echo "------------------------------------------------------------------------"
				echo
				python computeMeanAtQuantile.py ${var2D_ids} -p ${pr_ids} -q ${percentiles} \
				 -d ${dataroot} -c $compset -e $experiment -s $subset
			done
		done
	done
fi

## Compute mean profiles at locations of precipitation quantile
computeAllMeanProfilesAtQuantile=false
if [ "$computeAllMeanProfilesAtQuantile" == "true" ]
then
	for compset in `echo $compsets`
	do
		for experiment in `echo $experiments`
		do
			for subset in `echo $subsets`
			do
				echo 
				echo "------------------------------------------------------------------------"
				echo "     Computing mean profiles of 3D variables for:"
				echo "     - compset $compset"
				echo "     - experiment $experiment"
				echo "     - subset $subset"
				echo "     - target variables ${var3D_ids}"
				echo "     over location of percentiles $percentiles"
				echo "     for precipitation variables ${pr_ids}"
				echo "------------------------------------------------------------------------"
				echo
				python computeMeanProfileAtQuantile.py ${var3D_ids} -p ${pr_ids} -q ${percentiles} \
				 -d ${dataroot} -c $compset -e $experiment -s $subset
			done
		done
	done
fi

cd -