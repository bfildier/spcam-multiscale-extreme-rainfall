#!/bin/bash

## action can be either "preprocessing", "plotting"
action=$1
echo 
echo " -- Load environment for machine "$HOSTNAME" for "$1"..."

if [[ "$HOSTNAME" == "jollyjumper" ]]
then

	echo "Environment already loaded"

elif [[ "$HOSTNAME" == "edison"* ]]
then
	echo "- Loading modules..."
	
	if [ "$action" == "analysis" ]
	then 

		module load python/2.7-anaconda

	elif [ "$action" == "preprocessing" ]
	then

		module load python/2.7-anaconda
		module load nco

	fi
fi