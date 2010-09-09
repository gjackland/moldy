#!/bin/bash

## Run steps version 1.1	##
## by Martin Uhrin		##

# Turn on debugging
#set -x

## File/directory settings ##

readonly tools_reldir="../../../tools/"		# Path to directory where all the tools binaries are stored
readonly moldy_exe="../../../build/moldy.noOMP.FeC"			# The moldy binary to use for running the sim, must be compiled for right potential
readonly preproc_file="preprocess.sh"		# Preprocess file that can be executed to perform step specific actions
readonly input_dir="../system.in"			# Directory with system input files
readonly input_files=("../system.in/15^3.in") 	#$input_dir/*		# Input files to use (for all files in input_dir use $input_dir/*)
readonly atoms_token="nm"				# Token used in params.in to designate number of atoms
readonly timing_out="timing.out"			# The file where timing data will be stored
readonly runtime_string="Run time:"		# The string that preceeds the execution time as printed by MOLDY

## Run settings ##
readonly num_threads=(1)	# OMP_NUM_THREADS numbers to use
readonly num_repeats=3		# Number of times to repeat runs for each set of variables

readonly my_wd="`pwd`"

# Find the absolute path to the bin dir so that preprocess scripts can also use it
bin_dir="`pwd`/$tools_reldir"
if [ -a $bin_dir ]; then
	cd "$bin_dir"
	bin_dir="`pwd`"
	export BIN_DIR="$bin_dir"
	cd "$my_wd"
else
	echo "Error: directory with tools binaries ($bin_dir) does not exist"
	exit 1
fi

# Declare variables we will use
declare -i num_atoms=0
declare -i line_no=0

echo "Binaries directory is $bin_dir"

starttime="`date +%H:%M:%S-%d/%m/%y`"
echo -e "\n## Starting new run at $starttime ##" >> $timing_out
echo -e "# input_file\tnum_atoms\tnum_threads\ttimestamp\truntime\t...(repeated runs)" >> $timing_out

for input in $input_files
do
	echo "Processing input $input"
	
	# Let's make sure the number of atoms in system.in and params.in agree
	java -classpath "$tools_reldir" FixParams $input params.in

	# Get the number of atoms from the system input file
	num_atoms=`grep "$atoms_token" params.in | grep -m 1 -Eo [[:digit:]]+`

	if [ -a $preproc_file ]; then
		echo -e "Preprocessing"
		./$preproc_file
		echo -e "Finished preprocessing"
	fi
	
	for i in ${num_threads[@]}
	do
		echo -e "\tProcessing OMP_NUM_THREADS=$i"
		export OMP_NUM_THREADS=$i
		# Print the input file used, number of atoms and number of threads to the timing file
		echo -ne "$input\t$num_atoms\t$i" >> $timing_out

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
