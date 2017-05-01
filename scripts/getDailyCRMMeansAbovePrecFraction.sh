#!/bin/bash

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Load modules if necessary
source $scriptdir/loadPlatformSpecific.sh preprocessing

#---------------------------#
#   Experiment variables    #
#---------------------------#

compset="FSPCAMm_AMIP"
experiment="abrupt4xCO2"
case="bf_"${compset}"_"$experiment

if [[ "$HOSTNAME" == "jollyjumper" ]]
then
	inputdir="/Users/bfildier/Data/simulations/"${case}
elif [[ "$HOSTNAME" == "edison"* ]]
then	
	inputdir=${scriptdir%\/*}"/archive/"${case}"/atm/hist"
fi

#-----------------------------------#
#   Time boundaries and functions   #
#-----------------------------------#

lowerTimeBound="1850-05-01-03600"
upperTimeBound="1850-05-02-00000"

## Get start and end date
lowerDate=${lowerTimeBound%-*}
upperDate=${upperTimeBound%-*}

## Use start and end times for dividing all days
lowerTime=${lowerTimeBound##*-}
upperTime=${upperTimeBound##*-}

## Define function to increment date in correct format
cat << EOF > incrementDate
#!/usr/bin/python
import datetime as dt
import sys
newDate = dt.datetime.strptime(sys.argv[1],'%Y-%m-%d')
newDate += dt.timedelta(days=1)
print newDate.isoformat().split("T")[0]
EOF
chmod +x incrementDate

#-------------------------------------------------------#
#    Compute means for a sequence of fraction values    #
#-------------------------------------------------------#

# varids="CRM_T CRM_W CRM_QC CRM_QPC CRM_PREC"
varids="CRM_PREC"
# fractions="90 75 50 25 10"
fractions="100"
for varid in `echo $varids`
do 
	for fraction in `echo $fractions`
	do
		echo 
		echo "----------------------------------------------------------------------------------"
		echo "     Computing convective-scale means for "$varid" and with fraction "$fraction"% "
		echo "     for $compset and $experiment"
		echo "----------------------------------------------------------------------------------"
		startDate=$lowerDate
		while [ "$startDate" != "$upperDate" ]
		do
			endDate=`./incrementDate $startDate`
			# std message
			echo
			echo "$(tput bold) Extract values for ${startDate}-${lowerTime}"\
			"--- ${endDate}-${upperTime}$(tput sgr0)"
			echo

			python getDailyCRMMeansAbovePrecFraction.py $varid $fraction $experiment \
			           $case ${startDate}-${lowerTime} ${endDate}-${upperTime} $inputdir

			# Increment dates
			startDate=$endDate
		done
	done
done

exit 0