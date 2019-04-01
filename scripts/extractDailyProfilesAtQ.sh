#!/bin/bash

## Here execute python script extractDailyProfilesAtQ.py for:
## - all dates
## - a given value of:
##    . a 3D variable varid
##    . pr_id
##    . q_id
##    . experiment
##    . compset
##    . subset

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Load modules if necessary
source $scriptdir/loadPlatformSpecific.sh preprocessing

#---------------------------#
#   Experiment variables    #
#---------------------------#

compset="FSPCAMm_AMIP"
experiment="piControl"
case="bf_"${compset}"_"$experiment
subset="tropics"
varid="CRM_W_O50"
pr_id="PRECT"
q_id="98.7411"

if [[ "$HOSTNAME" == "jollyjumper" ]]
then
	dataroot="/Users/bfildier/Data"
elif [[ "$HOSTNAME" == "edison"* ]] || [[ "$HOSTNAME" == "cori"* ]]
then	
	dataroot=${scriptdir%\/*}
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

#--------------------------------------#
#    Extract indices of percentiles----#
#--------------------------------------#

startDate=$lowerDate

while [ "$startDate" != "$upperDate" ]
do
	endDate=`./incrementDate $startDate`
	# std message
	echo "$(tput bold) ${startDate}-${lowerTime} --- ${endDate}-${upperTime} $(tput sgr0)"

	python extractDailyProfilesAtQ.py -v ${varid} -p ${pr_id} -q ${q_id} -d ${dataroot} \
		-e ${experiment} -s ${subset} -c ${compset} -dt ${startDate}-${lowerTime} \
		${endDate}-${upperTime}

	# Increment dates
	startDate=$endDate
done

echo

exit 0
