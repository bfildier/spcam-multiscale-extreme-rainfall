#!/bin/bash

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $scriptdir

#-----------------------------#
#   Compute all statistics    #
#-----------------------------#

if [[ "$HOSTNAME" == "jollyjumper" ]]
then
	dataroot="/Users/bfildier/Data"
elif [[ "$HOSTNAME" == "edison"* ]]
then	
	dataroot=${scriptdir%\/*}
fi

compsets="FSPCAMm_AMIP FAMIPC5"
experiments="piControl abrupt4xCO2"
varids="PRECT PRECSH PRECL PRECC"
lowerDate="1850-05-01-00000"
upperDate="1850-06-01-00000"

for compset in `echo $compsets`
do
	for experiment in `echo $experiments`
	do
		for varid in `echo $varids`
		do

			echo 
			echo "------------------------------------------------------------------------"
			echo "     Extracting hourly $varid between $lowerDate and $upperDate for:"
			echo "     - compset $compset"
			echo "     - experiment $experiment"
			echo "------------------------------------------------------------------------"
			echo

			sed -i'' -e "s/^experiment=.*/experiment=\"${experiment}\"/" extractVar.sh
			sed -i'' -e "s/^compset=.*/compset=\"${compset}\"/" extractVar.sh
			sed -i'' -e "s/^lowerDate=.*/lowerDate=\"${lowerDate}\"/" extractVar.sh
			sed -i'' -e "s/^upperDate=.*/upperDate=\"${upperDate}\"/" extractVar.sh

			./extractVar.sh $varid

		done
	done
done



exit 0