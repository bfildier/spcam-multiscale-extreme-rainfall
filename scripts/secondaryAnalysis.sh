#!/bin/bash

## Here execute analysis that both depends on previous analysis results (such as
## precipitation statistics) and requires looking at the initial raw data (e.g.
## CRM-level variables at given precipitation percentiles)

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $scriptdir

# Load modules if necessary
source $scriptdir/loadPlatformSpecific.sh analysis

## Variables
compsets="FSPCAMm_AMIP FAMIPC5"
compsets_crm="FSPCAMm_AMIP"
experiments="piControl abrupt4xCO2"
subsets="tropics"
# pr_ids="PRECT CRM_PREC_I90 CRM_PREC_I75 CRM_PREC_I50 CRM_PREC_I25 CRM_PREC_I10"
pr_ids="PRECT CRM_PREC_I50"
area_ids="I50"
# percentiles="90.0      92.0567  93.6904  94.9881  96.0189  96.8377  97.4881  98.0047"\
# "  98.4151  98.7411  99.0      99.2057  99.369   99.4988  99.6019  99.6838"\
# "  99.7488  99.8005  99.8415  99.8741  99.9     99.9206  99.9369  99.9499"\
# "  99.9602  99.9684  99.9749  99.98    99.9842  99.9874  99.99    99.9921"\
# "  99.9937  99.995   99.996   99.9968  99.9975  99.998   99.9984  99.9987"
# percentiles="99.9"
percentiles="90.0      92.0567  93.6904  94.9881  96.0189  96.8377  97.4881  98.0047"\
"  98.4151  98.7411  99.0      99.2057  99.369   99.4988  99.6019  99.6838"\
"  99.7488  99.8005  99.8415  99.8741  99.9206  99.9369  99.9499"\
"  99.9602  99.9684  99.9749  99.98    99.9842  99.9874  99.99    99.9921"\
"  99.9937  99.995   99.996   99.9968  99.9975  99.998   99.9984  99.9987"
# var3D_ids="OMEGA_CRM_W_I90 RHO_CRM_T_I90 CRM_OMEGA_WT_I90 QVSATENV_CRM_T_I90"\
# " OMEGA_CRM_W_I75 RHO_CRM_T_I75 CRM_OMEGA_WT_I75 QVSATENV_CRM_T_I75"\
# " OMEGA_CRM_W_I25 RHO_CRM_T_I25 CRM_OMEGA_WT_I25 QVSATENV_CRM_T_I25"\
# " OMEGA_CRM_W_I10 RHO_CRM_T_I10 CRM_OMEGA_WT_I10 QVSATENV_CRM_T_I10"
var3D_ids="CRM_W_O50"

## Time bounds
lowerTimeBound="1850-05-01-03600"
if [[ "$HOSTNAME" == "jollyjumper" ]]
then
	upperTimeBound="1850-05-02-00000"
elif [[ "$HOSTNAME" == "edison"* ]] || [[ "$HOSTNAME" == "cori"* ]]
then
	upperTimeBound="1851-05-01-00000"
fi

###--- Main code ---###

## Extract GCM indices of given percentiles
extractGCMindicesOfPercentiles=false
if [ "$extractGCMindicesOfPercentiles" == "true" ]; then

	sed -i'' -e "s/^lowerTimeBound=.*/lowerTimeBound=\"${lowerTimeBound}\"/" getDailyGCMIndOfPercentile.sh
	sed -i'' -e "s/^upperTimeBound=.*/upperTimeBound=\"${upperTimeBound}\"/" getDailyGCMIndOfPercentile.sh
			
	for compset in `echo $compsets`; do
		for experiment in `echo $experiments`; do
			for subset in `echo $subsets`; do
				for pr_id in `echo ${pr_ids}`; do
					for percentile in `echo $percentiles`; do

						echo 
						echo "------------------------------------------------------------------------"
						echo "     Extracting GCM indices for:"
						echo "     - compset $compset"
						echo "     - experiment $experiment"
						echo "     - subset $subset"
						echo "     - precipitation variables ${pr_id}"
						echo "     - percentile $percentile"
						echo "------------------------------------------------------------------------"
						echo

						sed -i'' -e "s/^compset=.*/compset=\"${compset}\"/" getDailyGCMIndOfPercentile.sh
						sed -i'' -e "s/^experiment=.*/experiment=\"${experiment}\"/" getDailyGCMIndOfPercentile.sh
						sed -i'' -e "s/^subset=.*/subset=\"${subset}\"/" getDailyGCMIndOfPercentile.sh
						sed -i'' -e "s/^pr_id=.*/pr_id=\"${pr_id}\"/" getDailyGCMIndOfPercentile.sh
						sed -i'' -e "s/^percentile=.*/percentile=\"${percentile}\"/" getDailyGCMIndOfPercentile.sh
						./getDailyGCMIndOfPercentile.sh

					done
				done
			done
		done
	done
fi

## Extract GCM indices of all extreme percentiles simultaneously
extractGCMindicesAllExtremes=false
if [ "$extractGCMindicesAllExtremes" == "true" ]; then

	sed -i'' -e "s/^lowerTimeBound=.*/lowerTimeBound=\"${lowerTimeBound}\"/" getDailyGCMIndAllExtremes.sh
	sed -i'' -e "s/^upperTimeBound=.*/upperTimeBound=\"${upperTimeBound}\"/" getDailyGCMIndAllExtremes.sh
			
	for compset in `echo $compsets`; do
		for experiment in `echo $experiments`; do
			for subset in `echo $subsets`; do
				for pr_id in `echo ${pr_ids}`; do

					echo 
					echo "------------------------------------------------------------------------"
					echo "     Extracting GCM indices for:"
					echo "     - compset $compset"
					echo "     - experiment $experiment"
					echo "     - subset $subset"
					echo "     - precipitation variables ${pr_id}"
					echo "------------------------------------------------------------------------"
					echo

					sed -i'' -e "s/^compset=.*/compset=\"${compset}\"/" getDailyGCMIndAllExtremes.sh
					sed -i'' -e "s/^experiment=.*/experiment=\"${experiment}\"/" getDailyGCMIndAllExtremes.sh
					sed -i'' -e "s/^subset=.*/subset=\"${subset}\"/" getDailyGCMIndAllExtremes.sh
					sed -i'' -e "s/^pr_id=.*/pr_id=\"${pr_id}\"/" getDailyGCMIndAllExtremes.sh
					sed -i'' -e "s/^percentile=.*/percentile=\"${percentile}\"/" getDailyGCMIndAllExtremes.sh
					./getDailyGCMIndAllExtremes.sh

				done
			done
		done
	done
fi

## Extract vertically averaged CRM mass flux of given percentiles
extractCRMMFatQ=false
if [ "$extractCRMMFatQ" == "true" ]; then

	sed -i'' -e "s/^lowerTimeBound=.*/lowerTimeBound=\"${lowerTimeBound}\"/" extractCRMMFatQ.sh
	sed -i'' -e "s/^upperTimeBound=.*/upperTimeBound=\"${upperTimeBound}\"/" extractCRMMFatQ.sh
			
	for compset in `echo ${compsets_crm}`; do
		for experiment in `echo $experiments`; do
			for subset in `echo $subsets`; do
				for pr_id in `echo ${pr_ids}`; do
					for area_id in `echo ${area_ids}`; do
						for percentile in `echo ${percentiles}`; do

							echo 
							echo "------------------------------------------------------------------------"
							echo "     Extracting CRM mass flux for:"
							echo "     - compset $compset"
							echo "     - experiment $experiment"
							echo "     - subset $subset"
							echo "     - precipitation variables ${pr_id}"
							echo "     - percentile $percentile"
							echo "------------------------------------------------------------------------"
							echo

							sed -i'' -e "s/^compset=.*/compset=\"${compset}\"/" extractCRMMFatQ.sh
							sed -i'' -e "s/^experiment=.*/experiment=\"${experiment}\"/" extractCRMMFatQ.sh
							sed -i'' -e "s/^subset=.*/subset=\"${subset}\"/" extractCRMMFatQ.sh
							sed -i'' -e "s/^pr_id=.*/pr_id=\"${pr_id}\"/" extractCRMMFatQ.sh
							sed -i'' -e "s/^area_id=.*/area_id=\"${area_id}\"/" extractCRMMFatQ.sh
							sed -i'' -e "s/^q_id=.*/q_id=\"${percentile}\"/" extractCRMMFatQ.sh
							./extractCRMMFatQ.sh

						done
					done
				done
			done
		done
	done
fi

## Extract vertically averaged CRM mass flux of given percentiles in the dry area
extractCRMMFDryAtQ=false
if [ "$extractCRMMFDryAtQ" == "true" ]; then

	sed -i'' -e "s/^lowerTimeBound=.*/lowerTimeBound=\"${lowerTimeBound}\"/" extractCRMMFDryAtQ.sh
	sed -i'' -e "s/^upperTimeBound=.*/upperTimeBound=\"${upperTimeBound}\"/" extractCRMMFDryAtQ.sh
			
	for compset in `echo ${compsets_crm}`; do
		for experiment in `echo $experiments`; do
			for subset in `echo $subsets`; do
				for pr_id in `echo ${pr_ids}`; do
					for area_id in `echo ${area_ids}`; do
						for percentile in `echo ${percentiles}`; do

							echo 
							echo "------------------------------------------------------------------------"
							echo "     Extracting CRM mass flux for:"
							echo "     - compset $compset"
							echo "     - experiment $experiment"
							echo "     - subset $subset"
							echo "     - precipitation variables ${pr_id}"
							echo "     - percentile $percentile"
							echo "------------------------------------------------------------------------"
							echo

							sed -i'' -e "s/^compset=.*/compset=\"${compset}\"/" extractCRMMFDryAtQ.sh
							sed -i'' -e "s/^experiment=.*/experiment=\"${experiment}\"/" extractCRMMFDryAtQ.sh
							sed -i'' -e "s/^subset=.*/subset=\"${subset}\"/" extractCRMMFDryAtQ.sh
							sed -i'' -e "s/^pr_id=.*/pr_id=\"${pr_id}\"/" extractCRMMFDryAtQ.sh
							sed -i'' -e "s/^area_id=.*/area_id=\"${area_id}\"/" extractCRMMFDryAtQ.sh
							sed -i'' -e "s/^q_id=.*/q_id=\"${percentile}\"/" extractCRMMFDryAtQ.sh
							./extractCRMMFDryAtQ.sh

						done
					done
				done
			done
		done
	done
fi


## Extract vertical convective-scale profiles of W of given percentiles
extractCRMWatQ=false
if [ "$extractCRMWatQ" == "true" ]; then

	sed -i'' -e "s/^lowerTimeBound=.*/lowerTimeBound=\"${lowerTimeBound}\"/" extractCRMWatQ.sh
	sed -i'' -e "s/^upperTimeBound=.*/upperTimeBound=\"${upperTimeBound}\"/" extractCRMWatQ.sh
			
	for compset in `echo ${compsets_crm}`; do
		for experiment in `echo $experiments`; do
			for subset in `echo $subsets`; do
				for pr_id in `echo ${pr_ids}`; do
					for area_id in `echo ${area_ids}`; do
						for percentile in `echo ${percentiles}`; do

							echo 
							echo "------------------------------------------------------------------------"
							echo "     Extracting CRM mass flux for:"
							echo "     - compset $compset"
							echo "     - experiment $experiment"
							echo "     - subset $subset"
							echo "     - precipitation variables ${pr_id}"
							echo "     - percentile $percentile"
							echo "------------------------------------------------------------------------"
							echo

							sed -i'' -e "s/^compset=.*/compset=\"${compset}\"/" extractCRMWatQ.sh
							sed -i'' -e "s/^experiment=.*/experiment=\"${experiment}\"/" extractCRMWatQ.sh
							sed -i'' -e "s/^subset=.*/subset=\"${subset}\"/" extractCRMWatQ.sh
							sed -i'' -e "s/^pr_id=.*/pr_id=\"${pr_id}\"/" extractCRMWatQ.sh
							sed -i'' -e "s/^area_id=.*/area_id=\"${area_id}\"/" extractCRMWatQ.sh
							sed -i'' -e "s/^q_id=.*/q_id=\"${percentile}\"/" extractCRMWatQ.sh
							./extractCRMWatQ.sh

						done
					done
				done
			done
		done
	done
fi


## Extract vertical profiles of varid at locations of given percentiles
extractDailyProfilesAtQ=true
if [ "$extractDailyProfilesAtQ" == "true" ]; then

	sed -i'' -e "s/^lowerTimeBound=.*/lowerTimeBound=\"${lowerTimeBound}\"/" extractDailyProfilesAtQ.sh
	sed -i'' -e "s/^upperTimeBound=.*/upperTimeBound=\"${upperTimeBound}\"/" extractDailyProfilesAtQ.sh
			
	for compset in `echo ${compsets_crm}`; do
		for experiment in `echo $experiments`; do
			for subset in `echo $subsets`; do
				for varid in `echo ${var3D_ids}`; do
					for pr_id in `echo ${pr_ids}`; do
						for percentile in `echo ${percentiles}`; do

							echo 
							echo "------------------------------------------------------------------------"
							echo "     Extracting CRM mass flux for:"
							echo "     - compset $compset"
							echo "     - experiment $experiment"
							echo "     - subset $subset"
							echo "     - varid $varid"
							echo "     - precipitation variables ${pr_id}"
							echo "     - percentile $percentile"
							echo "------------------------------------------------------------------------"
							echo

							sed -i'' -e "s/^compset=.*/compset=\"${compset}\"/" extractDailyProfilesAtQ.sh
							sed -i'' -e "s/^experiment=.*/experiment=\"${experiment}\"/" extractDailyProfilesAtQ.sh
							sed -i'' -e "s/^subset=.*/subset=\"${subset}\"/" extractDailyProfilesAtQ.sh
							sed -i'' -e "s/^pr_id=.*/pr_id=\"${pr_id}\"/" extractDailyProfilesAtQ.sh
							sed -i'' -e "s/^varid=.*/varid=\"${varid}\"/" extractDailyProfilesAtQ.sh
							sed -i'' -e "s/^q_id=.*/q_id=\"${percentile}\"/" extractDailyProfilesAtQ.sh
							./extractDailyProfilesAtQ.sh

						done
					done
				done
			done
		done
	done
fi







			
