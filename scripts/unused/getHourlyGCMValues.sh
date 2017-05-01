#!/bin/bash

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Load modules if necessary
source $scriptdir/loadPlatformSpecific.sh plotting

#---------------------------#
#   Experiment variables    #
#---------------------------#

compset="FAMIPC5"
experiment="abrupt4xCO2"
case="bf_"${compset}"_"$experiment

if [[ "$HOSTNAME" == "jollyjumper" ]]
then
	inputdir="/Users/bfildier/Data/simulations/"${case}
elif [[ "$HOSTNAME" == "edison"* ]]
then	
	inputdir=${scriptdir%\/*}"/archive/"${case}
fi

#----------------------#
#   Time boundaries    #
#----------------------#

lowerTimeBound="1850-05-01-03600"
upperTimeBound="1850-05-02-00000"

#-----------------------------------------#
#    Get values for a set of variables    #
#-----------------------------------------#

varids="CAPE CBMF CIN EVAPPREC FU FV ICEFRAC OMEGA PCLDBOT PCLDTOP PDELDRY"\
" PRECC PRECL PRECSH PRECT PRECTMX PS Q RELHUM SPMC SPMCUP SST T TMQ TS"\
" U V cin_Cu cinlcl_Cu umf_Cu wu_Cu"
# varids="PDELDRY PRECC TS U"
for varid in `echo $varids`
do
	echo
	echo "------------------------------------------------------------------"
	echo "     Getting hourly $varid values for $compset and $experiment"
	echo "------------------------------------------------------------------"
	python getHourlyGCMValues.py $varid $experiment $case $lowerTimeBound $upperTimeBound $inputdir
done

exit 0
