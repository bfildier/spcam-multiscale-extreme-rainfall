#!/bin/bash

## Here execute python script getDailyGCMIngAllExtremes.py for:
## - all dates
## - all percentiles
## - a given value of:
##    . pr_id
##    . experiment
##    . compset
##    . subset

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Load modules if necessary
source $scriptdir/loadPlatformSpecific.sh preprocessing

#---------------------------#
#   Experiment variables    #
#---------------------------#

compset="FAMIPC5"
experiment="abrupt4xCO2"
case="bf_"${compset}"_"$experiment
subset="tropics"
pr_id="CRM_PREC_I10"

if [[ "$HOSTNAME" == "jollyjumper" ]]
then
	dataroot="/Users/bfildier/Data"
elif [[ "$HOSTNAME" == "edison"* ]]
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

# echo 
# echo "$(tput bold)     Extracting GCM indices of "$Q"th percentile of "$pr_id" in subset "$subset"   $(tput sgr0)"
# echo "$(tput bold)     for $compset and $experiment$(tput sgr0)"
# echo 

while [ "$startDate" != "$upperDate" ]
do
	endDate=`./incrementDate $startDate`
	# std message
	echo "$(tput bold) ${startDate} --- ${endDate} $(tput sgr0)"

	python getDailyGCMIndAllExtremes.py -p ${pr_id} -d $dataroot -e $experiment \
	 -s $subset -c $compset -dt ${startDate} ${endDate}

	# Increment dates
	startDate=$endDate
done

echo

exit 0


