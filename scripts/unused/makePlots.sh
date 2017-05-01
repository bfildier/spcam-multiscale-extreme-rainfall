#!/bin/bash

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Load modules if necessary
source $scriptdir/loadPlatformSpecific.sh plotting

#---------------------#
#   Make all plots    #
#---------------------#

if [[ "$HOSTNAME" == "jollyjumper" ]]
then
	dataroot="/Users/bfildier/Data"
elif [[ "$HOSTNAME" == "edison"* ]]
then	
	dataroot=${scriptdir%\/*}
fi

##-- Plot fractional change delta P vs. extreme percentile Q --##

## Compare parameterized and superparameterized runs
fracChangeSPCAMvsCAM=false

## Plot one variable at a time
fracChange=false
if [ "$fracChange" == "true" ]
then
	# pr_varids="PRECT PRECC CRM_PREC_I90 CRM_PREC_I75 CRM_PREC_I50 CRM_PREC_I25 CRM_PREC_I10"
	pr_varids="PRECT PRECC CRM_PREC_I90 CRM_PREC_I75 CRM_PREC_I50"
	for pr_varid in `echo "${pr_varids}"`
	do
		echo
		echo " -- Plot fractional change in ${pr_varid} percentile P_q vs. q"
		python plotFractionalPvsExtremeQ.py ${pr_varid} ${dataroot}
	done
fi

# All on same graph
fracChangeAll=true
if [ "$fracChangeAll" == "true" ]
then
	# pr_varids="PRECT PRECC CRM_PREC_I90 CRM_PREC_I75 CRM_PREC_I50 CRM_PREC_I25 CRM_PREC_I10"
	pr_varids="PRECT CRM_PREC_I90 CRM_PREC_I75 CRM_PREC_I50 CRM_PREC_I25 CRM_PREC_I10"
	echo
	echo " -- Plot fractional change in ${pr_varids} percentile P_q vs. q"
	python plotFractionalPvsExtremeQ.py ${pr_varids} -d ${dataroot}
fi

##-- Plot mean vertical profiles of T, q, H, rho, omega and w --##
##-- for the whole Tropics or for extreme P quantile          --##

verticalProfiles=false
if [ "$verticalProfiles" == "true" ]
then
	pr_varids="PRECT PRECC CRM_PREC_I90 CRM_PREC_I75 CRM_PREC_I50 CRM_PREC_I25 CRM_PREC_I10"
	quantiles="90 95 99 99.9 99.99"
	for pr_varid in `echo ${pr_varids}`
	do
		for q in `echo $quantiles`
		do
			echo 
			echo " -- Plot vertical profiles for tropical mean and for "$q" extreme "${pr_varid}
			python plotProfilesForExtremeP.py ${pr_varid} $q ${dataroot}
		done
	done
fi