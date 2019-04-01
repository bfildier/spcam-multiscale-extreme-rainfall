#!/bin/bash

############################################################################
#                                                                          #
#   From CESM history files, extract in one single file a given variable   #
#   given in argument and store it in a specific folder.                   #
#                                                                          #
#                                               Ben Fildier, June 2017     #
#                                                                          #
############################################################################

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $scriptdir

# Load modules if necessary
source $scriptdir/loadPlatformSpecific.sh preprocessing

varid=$1

compset="FSPCAMm_AMIP"
experiment="piControl"
# experiment="abrupt4xCO2"
case="bf_"${compset}"_"$experiment

if [[ "$HOSTNAME" == "jollyjumper" ]]
then
	inputdir="/Users/bfildier/Data/simulations/"${case}
	outputdir=$inputdir/"extracted"
elif [[ "$HOSTNAME" == "edison"* ]]
then	
	inputdir=${scriptdir%\/*}"/archive/"${case}"/atm/hist"
	outputdir=${scriptdir%\/*}"/outputs/"${case}
fi

#----------------------#
#   Space boundaries   #
#----------------------#

getRegion='true'
latmin=7.0
latmax=9.0
lonmin=166.5
lonmax=169.0

#-----------------------------------#
#   Time boundaries and functions   #
#-----------------------------------#

lowerDate="1850-05-01-00000"
upperDate="1850-05-02-00000"

## Define function to increment date in correct format
cat << EOF > incrementDate
#!/usr/bin/python
import datetime as dt
import sys,string
newDateTime = dt.datetime.strptime(sys.argv[1][:10],'%Y-%m-%d')
delta = dt.timedelta(seconds=int(sys.argv[1].split('-')[-1]))
newDateTime += delta
newDateTime += dt.timedelta(hours=1)
newTime = dt.timedelta(hours=newDateTime.hour,minutes=newDateTime.minute,seconds=newDateTime.second)
print newDateTime.isoformat().split("T")[0]+'-'+"%05d"%newTime.seconds
print 
EOF
chmod +x incrementDate


#-----------------------------------------------#
#    Get values for a sequences of variables    #
#-----------------------------------------------#

echo $inputdir
cd $inputdir

echo
echo "------------------------------------------------------------------"
echo "     Getting hourly $varid values for $compset and $experiment"
echo "------------------------------------------------------------------"
startDate=$lowerDate
# global startDate
while [ "$startDate" != "$upperDate" ]
do
	endDate=`$scriptdir/incrementDate $startDate`
	# std message
	# echo
	# echo "$(tput bold) Compute mass flux for ${startDate}"\
	# "--- ${endDate}$(tput sgr0)"
	# echo	

	echo ${startDate}_$endDate	

	if [ "$getRegion" == "true" ]
	then
		# lmin=`bc -l <<< "$latmin"`
		# echo $lmin
		ncks -d lat,$latmin,$latmax -d lon,$lonmin,$lonmax -v $varid,date,datesec,ndbase,nsbase,nbdate,nbsec *.h0.$endDate.nc temp.$endDate.nc
	else
		ncks -v $varid,date,datesec,ndbase,nsbase,nbdate,nbsec *.h0.$endDate.nc temp.$endDate.nc
	fi

	# Increment dates
	startDate=$endDate
done

if [ "$getRegion" == "true" ]
then
	filename=${varid}_1hr-lon${lonmin}-${lonmax}-lat${latmin}-${latmax}_${compset}_${experiment}_r1i1p1_$lowerDate-$upperDate.nc
else
	filename=${varid}_1hr_${compset}_${experiment}_r1i1p1_$lowerDate-$upperDate.nc
fi

# Remove pre-existing output file in current directory
if [ -f $filename ]
then 
	rm $filename
fi

# Concat and remove temporary files
ncrcat -v $varid,date,datesec,ndbase,nsbase,nbdate,nbsec temp.*.nc $filename
rm temp.*.nc

# Create output directory
if [ ! -d $outputdir ]
then 
	mkdir $outputdir
fi

# Remove pre-existing output file
if [ -f $outputdir/$filename ]
then 
	rm $outputdir/$filename
fi

# Move file to output directory
mv $filename $outputdir"/"

cd -

exit 0