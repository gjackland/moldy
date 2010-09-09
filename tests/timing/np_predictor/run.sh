#!/bin/bash

## Run version 1.0		##
## by Martin Uhrin		##

# Turn on debugging
#set -x

## File/directory settings ##
readonly run_script="../../../tools/run_steps.sh"
readonly input_dir="../system.in"			# Directory with system input files
readonly input_files=( "../system.in/25^3.in" ) 	#$input_dir/*		# Input files to use (for all files in input_dir use $input_dir/*)


## Export all the settings ##
export moldy_exe="../../../build/moldy.FeC"	# The moldy binary to use for running the sim, must be compiled for right potential
export x_input_files=${input_files[*]}
export num_repeats=3
export x_num_threads="16"	# OMP_NUM_THREADS numbers to use

## Finally run the script ##
./$run_script
