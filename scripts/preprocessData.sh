#!/bin/bash

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $scriptdir

# Load modules if necessary
source $scriptdir/loadPlatformSpecific.sh preprocessing

#-----------------------------------------#
#    Execute all preprocessing scripts    #
#-----------------------------------------#

#----------------------#
#    GCM operations    #
#----------------------#

lowerTimeBound="1850-05-01-03600"
upperTimeBound="1850-05-02-00000"

compsets="FSPCAMm_AMIP FAMIPC5"
experiments="piControl abrupt4xCO2"

## Extract hourly data to daily files
hourlyDataInDailyFiles=false
if [ "$hourlyDataInDailyFiles" == "true" ]
then
	for compset in `echo $compsets`;
	do
		for experiment in `echo $experiments`;
		do
			sed -i'' -e "s/^compset=.*/compset=\"${compset}\"/" getHourlyGCMValuesDailyFiles.sh
			sed -i'' -e "s/^experiment=.*/experiment=\"${experiment}\"/" getHourlyGCMValuesDailyFiles.sh
			sed -i'' -e "s/^lowerTimeBound=.*/lowerTimeBound=\"${lowerTimeBound}\"/" getHourlyGCMValuesDailyFiles.sh
			sed -i'' -e "s/^upperTimeBound=.*/upperTimeBound=\"${upperTimeBound}\"/" getHourlyGCMValuesDailyFiles.sh
			./getHourlyGCMValuesDailyFiles.sh | tee logs/processData_hourlyDataInDailyFiles.log
		done
	done
fi

# ## Extract hourly data to hourly files
# hourlyDataInHourlyFiles=true
# if [ "$hourlyDataInHourlyFiles" == "true" ]
# then
# 	for compset in `echo $compsets`;
# 	do
# 		for experiment in `echo $experiments`;
# 		do
# 			sed -i '' -e "s/^compset=.*/compset=\"${compset}\"/" getHourlyGCMValues.sh
# 			sed -i '' -e "s/^experiment=.*/experiment=\"${experiment}\"/" getHourlyGCMValues.sh
# 			./getHourlyGCMValues.sh | tee logs/processData_hourlyDataInHourlyFiles.log
# 		done
# 	done
# fi

## Hourly vertical mass flux to daily files
hourlyMassFluxToDailyFiles=false
if [ "$hourlyMassFluxToDailyFiles" == "true" ]
then
	for compset in `echo $compsets`;
	do
		for experiment in `echo $experiments`;
		do
			sed -i'' -e "s/^compset=.*/compset=\"${compset}\"/" getHourlyMassfluxDailyFiles.sh
			sed -i'' -e "s/^experiment=.*/experiment=\"${experiment}\"/" getHourlyMassfluxDailyFiles.sh
			sed -i'' -e "s/^lowerTimeBound=.*/lowerTimeBound=\"${lowerTimeBound}\"/" getHourlyMassfluxDailyFiles.sh
			sed -i'' -e "s/^upperTimeBound=.*/upperTimeBound=\"${upperTimeBound}\"/" getHourlyMassfluxDailyFiles.sh
			./getHourlyMassfluxDailyFiles.sh | tee logs/processData_hourlyMassFluxToDailyFiles.log
		done
	done
fi

## Extract daily means to daily files
dailyMeansToDailyFiles=false
if [ "$dailyMeansToDailyFiles" == "true" ]
then
	for compset in `echo $compsets`;
	do
		for experiment in `echo $experiments`;
		do
			sed -i'' -e "s/^compset=.*/compset=\"${compset}\"/" getDailyFromHourlyGCMValues.sh
			sed -i'' -e "s/^experiment=.*/experiment=\"${experiment}\"/" getDailyFromHourlyGCMValues.sh
			sed -i'' -e "s/^lowerTimeBound=.*/lowerTimeBound=\"${lowerTimeBound}\"/" getDailyFromHourlyGCMValues.sh
			sed -i'' -e "s/^upperTimeBound=.*/upperTimeBound=\"${upperTimeBound}\"/" getDailyFromHourlyGCMValues.sh
			./getDailyFromHourlyGCMValues.sh | tee logs/processData_dailyMeansToDailyFiles.log
		done
	done
fi


#----------------------#
#    CRM operations    #
#----------------------#

compsets="FSPCAMm_AMIP"
experiments="piControl abrupt4xCO2"
# experiments="abrupt4xCO2"

## Extract hourly indices of highest XX% local CRM_PREC values
dailyCRMIndicesToDailyFiles=false
if [ "$dailyCRMIndicesToDailyFiles" == "true" ]
then
	for compset in `echo $compsets`;
	do
		for experiment in `echo $experiments`;
		do
			sed -i'' -e "s/^compset=.*/compset=\"${compset}\"/" getDailyCRMIndAbovePrecFraction.sh
			sed -i'' -e "s/^experiment=.*/experiment=\"${experiment}\"/" getDailyCRMIndAbovePrecFraction.sh
			sed -i'' -e "s/^lowerTimeBound=.*/lowerTimeBound=\"${lowerTimeBound}\"/" getDailyCRMIndAbovePrecFraction.sh
			sed -i'' -e "s/^upperTimeBound=.*/upperTimeBound=\"${upperTimeBound}\"/" getDailyCRMIndAbovePrecFraction.sh
			./getDailyCRMIndAbovePrecFraction.sh | tee logs/processData_dailyCRMIndicesToDailyFiles.log
		done
	done
fi

## Compute means for indices of highest XX% local CRM_PREC values
dailyMeansAbovePrecFraction=false
if [ "$dailyMeansAbovePrecFraction" == "true" ]
then
	for compset in `echo $compsets`;
	do
		for experiment in `echo $experiments`;
		do
			sed -i'' -e "s/^compset=.*/compset=\"${compset}\"/" getDailyCRMMeansAbovePrecFraction.sh
			sed -i'' -e "s/^experiment=.*/experiment=\"${experiment}\"/" getDailyCRMMeansAbovePrecFraction.sh
			sed -i'' -e "s/^lowerTimeBound=.*/lowerTimeBound=\"${lowerTimeBound}\"/" getDailyCRMMeansAbovePrecFraction.sh
			sed -i'' -e "s/^upperTimeBound=.*/upperTimeBound=\"${upperTimeBound}\"/" getDailyCRMMeansAbovePrecFraction.sh
			./getDailyCRMMeansAbovePrecFraction.sh | tee logs/processData_dailyMeansAbovePrecFraction_${compset}_${experiment}.log
		done
	done
fi

## Compute mean updrafts and downdrafts for indices of highest XX% local CRM_PREC values
dailyMeanDraftsPrecFraction=true
if [ "$dailyMeanDraftsPrecFraction" == "true" ]
then
	for compset in `echo $compsets`;
	do
		for experiment in `echo $experiments`;
		do
			sed -i'' -e "s/^compset=.*/compset=\"${compset}\"/" getDailyCRMMeanDraftsPrecFraction.sh
			sed -i'' -e "s/^experiment=.*/experiment=\"${experiment}\"/" getDailyCRMMeanDraftsPrecFraction.sh
			sed -i'' -e "s/^lowerTimeBound=.*/lowerTimeBound=\"${lowerTimeBound}\"/" getDailyCRMMeanDraftsPrecFraction.sh
			sed -i'' -e "s/^upperTimeBound=.*/upperTimeBound=\"${upperTimeBound}\"/" getDailyCRMMeanDraftsPrecFraction.sh
			./getDailyCRMMeanDraftsPrecFraction.sh | tee logs/processData_dailyMeanDraftsPrecFraction_${compset}_${experiment}.log
		done
	done
fi

## Compute daily mean area and fraction of convective rainfall events (PRECT/CRM_PREC_IXX)
rainEventGeometry=false
if [ "$rainEventGeometry" == "true" ]
then
	for compset in `echo $compsets`;
	do
		for experiment in `echo $experiments`;
		do
			sed -i'' -e "s/^compset=.*/compset=\"${compset}\"/" getRainEventGeometry.sh
			sed -i'' -e "s/^experiment=.*/experiment=\"${experiment}\"/" getRainEventGeometry.sh
			./getRainEventGeometry.sh | tee logs/processData_rainEventGeometry_${compset}_${experiment}.log
		done
	done
fi


exit 0