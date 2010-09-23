#!/bin/bash

## Run steps version 1.31	##
## by Martin Uhrin		##
#
# This file take a bunch of parameter inputs taken from export commands
# and uses them to run the timing test.
#

# Turn on debugging
#set -x

## File/directory settings ##

readonly tools_reldir="."				# Path to directory where all the tools binaries are stored
readonly atoms_token="nm"				# Token used in params.in to designate number of atoms
readonly timing_out="timing.out"			# The file where timing data will be stored

## Default values ##
declare -ir d_num_repeats=3				# Number of times to repeat each run
declare -ir d_num_threads=(1 8)			# Numbers of threads to use
declare -r d_runtime_string="MDloop runtime:"	# The string that preceeds the execution time as printed by MOLDY

## Check that required parameters are set ##
# The moldy binary to use for running the sim, must be compiled for right potential
if [ -z "${moldy_exe}" ]; then
	echo "Please export moldy_exe to point to MOLDY binary"
	exit 1
fi

# Input files to use 
declare -ar input_files=(${x_input_files})
if [ -z "${input_files[*]}" ]; then
	echo "Please set x_input_files to be a string of input files"
	exit 1
fi	

## use default values if not supplied ##

declare -ar num_threads=( ${x_num_threads=${d_num_threads[*]}} )	# OMP_NUM_THREADS numbers to use
declare -r runtime_string=${RUNTIME_STRING=${d_runtime_string}}	# String used to get time from MOLDY output

# Number of times to repeat runs for each set of variables
if [ -z "${num_repeats}" ]; then
	declare -ir num_repeats=$d_num_repeats
fi

###################################################
## Start of script actually doing stuff 		##

starttime="`date +%H:%M:%S-%d/%m/%y`"
echo -e "\n## Starting new run at $starttime ##" >> $timing_out

echo "--SETTINGS--" >> $timing_out
echo -e "Executable:\t$moldy_exe" >> $timing_out
echo -e "Repeats:\t$num_repeats" >> $timing_out
echo -e "Input files:\t${input_files[*]}" >> $timing_out

readonly my_wd="`pwd`"			# Get the current directory
readonly my_dir="`dirname $0`"	# Get the directory where the script is relative to wd

# Find the absolute path to the bin dir so that preprocess scripts can also use it
tools_dir="$my_dir/$tools_reldir"
if [ -a $tools_dir ]; then
	cd "$tools_dir"
	tools_dir="`pwd`"
	cd "$my_wd"
else
	echo "Error: directory with tools binaries ($tools_dir) does not exist"
	exit 1
fi

# Declare variables we will use
declare -i num_atoms=0
declare -i line_no=0


echo -ne "# input_file\tnum_atoms\tnum_threads" >> $timing_out
if [ -n "${custom_param}" ]; then
	echo -ne "\t${custom_param}" >> $timing_out
fi
echo -e "\ttimestamp\truntime\t...(repeated runs)" >> $timing_out

for input in ${input_files[@]}
do
	echo "Processing input $input"
	
	# Let's make sure the number of atoms in system.in and params.in agree
	java -classpath "$tools_dir" FixParams $input params.in

	# Get the number of atoms from the system input file
	num_atoms=`grep "$atoms_token" params.in | grep -m 1 -Eo [[:digit:]]+`
	
	for i in ${num_threads[@]}
	do
		echo -e "\tProcessing OMP_NUM_THREADS=$i"
		export OMP_NUM_THREADS=$i

		# Print the input file used, number of atoms and number of threads to the timing file
		echo -ne "$input\t$num_atoms\t$i" >> $timing_out

		# If the user has specified a custom parameter value then print that as well
		if [ -n "${custom_param}" ]; then
			echo -ne "\t$custom_param_value" >> $timing_out
		fi

		# Define the name of the file to store output from the simulation run
		output_file="N=${num_atoms}P=${i}.out"
	
		# Loop over repeats
		j="0"
		while [ $j -lt $num_repeats ]
		do
			# Put a timestamp in the output file
			run_timestamp="`date +%H:%M:%S-%d/%m/%y`"
			echo "== Output from run at $run_timestamp ==" >> $output_file

			line_no=`wc -l < $output_file`

			# Run MOLDY and append stdout and stderr to file
			$moldy_exe >> $output_file 2>> $output_file

			# Now extract the run time from the output file and save it
			runtime=`sed -n ${line_no},'$'p $output_file | grep "$runtime_string" | sed "s/$runtime_string//"`
			echo -ne "\t$run_timestamp\t$runtime" >> $timing_out

			j=$[$j+1]
		done

		# Put a newline at the end of this timing data
		echo "" >> $timing_out
		echo -e "\tFinished OMP_NUM_THREADS=$i"
	done

	echo "Finished input $input"
done
