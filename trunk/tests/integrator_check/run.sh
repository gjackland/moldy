#!/bin/bash

## Run parameters script 	##
## version 1.0		##
## by Martin Uhrin		##

# Turn debugging on
#set -x

## File/directory settings ##

readonly moldy_exe="../../build/moldy.FeC"
readonly tools_reldir="../../tools"			# Path to directory where all the tools binaries are stored
readonly param_set_script="replace_token.sh"		# The name of the script in the tools directory that sets the parameter to its new value
readonly param_token="iverlet"				# Token of the parameter to be varied
readonly param_values=(0 1)					# Values that the parameter will take

# Num threads to use
export OMP_NUM_THREADS=8

if [ ! -f $tools_reldir/$param_set_script ]; then
	echo "Parameter set script $tools_reldir/$param_set_script does not exist"
	exit
fi

for param_value in ${param_values[@]}
do
	out_file=${param_token}=${param_value}.out

	# Set the parameter in the params.in file
	sh $tools_reldir/$param_set_script params.in $param_token $param_value

	# Now start the run
	./$moldy_exe >> $out_file

	# Move the files of interest
	mv output66.txt output66${param_token}=${param_value}.txt
	mv system.out system${param_token}=${param_value}.out
done

