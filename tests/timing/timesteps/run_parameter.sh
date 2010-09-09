#!/bin/bash

## Run parameters script 	##
## version 1.0			##
## by Martin Uhrin		##

# Turn on debug
#set -x

## File/directory settings ##

readonly tools_reldir="../../../tools"		# Path to directory where all the tools binaries are stored
readonly run_script="run_steps.sh"			# Script to actually execute the timing runs
readonly param_set_script="replace_token.sh"	# The name of the script in the tools directory that sets the parameter to its new value
readonly param_token="nsteps"				# Token of the parameter to be varied
readonly param_values=(2500 7500 50000)		# Values that the 
parameter will take
readonly timing_out="timing.out"			# The file where timing data will be stored

if [ ! -f $tools_reldir/$param_set_script ]; then
	echo "Parameter set script '$tools_reldir/$param_set_script' does not exist"
	exit
fi

echo "PATEMETERS RUN: $param_token will take on values: ${param_values[*]}"

for param_value in ${param_values[@]}
do
	# Set the parameter in the params.in file
	sh $tools_reldir/$param_set_script params.in $param_token $param_value

	# Write the value of the parameter to the timing file
	echo "$param_token = $param_value" >> $timing_out

	# Now execute the timing run
	./$run_script
done

