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
	dataroot="/Users/bfildier/Data"
	inputdir="/Users/bfildier/Data/simulations/"${case}
elif [[ "$HOSTNAME" == "edison"* ]] || [[ "$HOSTNAME" == "cori"* ]]
then
	dataroot=${scriptdir%\/*}
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
# varids="QVSATENV_CRM_T RHO_CRM_T OMEGA_CRM_WT"
varids="CRM_OMEGA CRM_OMEGA_WT"
# varids="OMEGA_CRM_W RHO_CRM_T CRM_OMEGA_WT QVSATENV_CRM_T"
# fractions="90 75 25 10"
# fractions="-1"
fractions="50"
for varid in `echo $varids`
do 
	for fraction in `echo $fractions`
	do
		echo 
		echo "--------------------------------------------------------------------------------------"
		echo "     Computing convective-scale values of "$varid" and with fraction "$fraction"% "
		echo "     for $compset and $experiment"
		echo "--------------------------------------------------------------------------------------"
		startDate=$lowerDate
		while [ "$startDate" != "$upperDate" ]
		do
			endDate=`./incrementDate $startDate`
			# std message
			echo
			echo "$(tput bold) Extract values for ${startDate}-${lowerTime}"\
			"--- ${endDate}-${upperTime}$(tput sgr0)"
			echo

			if [ "$varid" == "CRM_OMEGA" ]		# Pressure velocity profile computed from hourly CRM-level W, hourly GCM-level RHO and IXX indices
			then 
				python computeDailyCRMOmegaAbovePrecFraction.py -f $fraction -e \
					$experiment -c $compset -d $dataroot -dt ${startDate}-${lowerTime} \
					${endDate}-${upperTime}
			elif [ "$varid" == "CRM_OMEGA_WT" ]		# Pressure velocity profile computed from hourly CRM-level W and T, hourly GCM-level Q and PS and indices IXX
			then
				python computeDailyCRMOmegaWTAbovePrecFraction.py -f $fraction -e \
					$experiment -c $compset -d $dataroot -dt ${startDate}-${lowerTime} \
					${endDate}-${upperTime}
			elif [ "$varid" == "QVSATENV_CRM_T" ] 	# Saturation specific humidity profile computed from convective-scale T (CRM_T_IXX) and mesoscale (daily) PS 
			then
				python computeDailyQvstarFromTAbovePrecFraction.py -f $fraction -e \
					$experiment -c $compset -d $dataroot -dt ${startDate}-${lowerTime} \
					${endDate}-${upperTime}
			elif [ "$varid" == "RHO_CRM_T" ] 	# Density profile computed from convective-scale T (CRM_T_IXX) and mesoscale (daily) PS and Q
			then
				python computeDailyDensityFromGCMQPAndTAbobePrecFraction.py -f $fraction -e \
					$experiment -c $compset -d $dataroot -dt ${startDate}-${lowerTime} \
					${endDate}-${upperTime}
			elif [ "$varid" == "OMEGA_CRM_WT" ]		# Pressure velocity profile computed from convective-scale RHO_CRM_T_IXX and CRM_W_IXX
			then
				python computeDailyOmegaFromWandTAbovePrecFraction.py -f $fraction -e \
					$experiment -c $compset -d $dataroot -dt ${startDate}-${lowerTime} \
					${endDate}-${upperTime}
			elif [ "$varid" == "OMEGA_CRM_W" ]		# Pressure velocity profile computed from convective-scale CRM_W_IXX and mesoscale RHO_MESO (itself computed from mesoscale T, Q, PS)
			# This requires having run this script beforehand with "$varid" == "RHO_CRM_T" and "$fraction" == "-1"
			then
				python computeDailyOmegaFromWAbovePrecFraction.py -f $fraction -e \
					$experiment -c $compset -d $dataroot -dt ${startDate}-${lowerTime} \
					${endDate}-${upperTime}
			else
				python getDailyCRMMeansAbovePrecFraction.py $varid $fraction $experiment \
			    	       $case ${startDate}-${lowerTime} ${endDate}-${upperTime} $inputdir
			fi

			# Increment dates
			startDate=$endDate
		done
	done
done

exit 0
