#!/bin/bash

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Load modules if necessary
source $scriptdir/loadPlatformSpecific.sh preprocessing

#---------------------------#
#   Experiment variables    #
#---------------------------#

compset="FAMIPC5"
experiment="abrupt4xCO2"
case="bf_"${compset}"_"$experiment

if [[ "$HOSTNAME" == "jollyjumper" ]]
then
	inputdir="/Users/bfildier/Data/preprocessed/"${case}"/1hr"
elif [[ "$HOSTNAME" == "edison"* ]] || [[ "$HOSTNAME" == "cori"* ]]
then	
	inputdir=${scriptdir%\/*}"/preprocessed/"${case}"/1hr"
fi

#---------------------------------------------------------------#
#    Convert hourly to daily values for a sequence of varids    #
#---------------------------------------------------------------#

# varids="CAPE CBMF CIN EVAPPREC FU FV ICEFRAC OMEGA PCLDBOT PCLDTOP PDELDRY"\
# " PRECC PRECL PRECSH PRECT PRECTMX PS Q RELHUM SPMC SPMCUP SST T TMQ TS"\
# " U V cin_Cu cinlcl_Cu"
varids="RHO"
for varid in `echo $varids`
do
	echo
	echo "------------------------------------------------------------"
	echo "------ Converting hourly to daily values for "$varid" ------"
	echo "------------------------------------------------------------"
	echo 
	for file in `ls ${inputdir}/${varid}_*`
	do
		dailyfile=`echo $file | sed -e 's/1hr/day/g' -e 's/0100-/-/' -e 's/0000.nc/.nc/'`
		echo ${file##*\/}
		echo "  -->  "${dailyfile##*\/}
		ncra -O -y avg $file $dailyfile
	done
done

exit 0
