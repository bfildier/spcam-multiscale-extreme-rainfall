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
percentiles="99.9"

## Time bounds
lowerTimeBound="1850-05-01-03600"
upperTimeBound="1851-05-01-00000"

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
extractCRMMFDryAtQ=true
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












			
