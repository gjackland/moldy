#!/bin/bash

## Correctness test version 1.1	##
## by Martin Uhrin			##

# Turn on debugging
#set -x

## File/directory settings ##

readonly tools_reldir="../../tools"			# Path to directory where all the tools binaries are stored
readonly moldy_exe="../../build/moldy.FeC"		# The moldy binary to use for running the sim, must be compiled for right potential
readonly reference_dir="../reference"			# Directory containing directories of reference results
readonly test_dirs=(`find $reference_dir/* -maxdepth 1 -type d | sed "s_${reference_dir}/__"`)	# Get all the test directories
readonly system_in_token="file_system"
readonly output_file="run.out"
readonly params_file="params.in"

readonly my_wd="`pwd`"

# Find the absolute path to the bin dir so that preprocess scripts can also use it
tools_dir="`pwd`/$tools_reldir"
if [ -a $bin_dir ]; then
	cd "$tools_dir"
	tools_dir="`pwd`"
	cd "$my_wd"
else
	echo "Error: directory with tools binaries ($bin_dir) does not exist"
	exit 1
fi

## Run settings ##
declare -ir sig_digits=11	# Results have to agree to this many significant digits
readonly num_threads=(8)	# OMP_NUM_THREADS numbers to use

# The command to comare two files
readonly compare="java -classpath $tools_dir:$tools_dir/args4j.jar Compare -d $sig_digits"

# Declare variables we will use

echo ${test_dirs[*]}

for test_dir in ${test_dirs[@]}
do

	echo "Processing reference test $test_dir"
	echo "$reference_dir/${test_dir}${params_file}"

	if [ ! -f "${reference_dir}/${test_dir}/${params_file}" ]; then
		echo "Skipping test $test, $params_file does not exist"
		continue
	fi

	# Get the system.in filename: first extract it and then strip off any quotes at start and end
	system_in=`sed -n "s/[[:blank:]]*${system_in_token}[[:blank:]]*=[[:blank:]]*\([^[:blank:]]*\).*/\1/p" \
			$reference_dir/${test_dir}/${params_file} | sed 's/"\(.*\)"/\1/'`

	# Now make sure the system input file exists
	if [ ! -f "${reference_dir}/${test_dir}/${system_in}" ]; then
		echo "Skipping test $test, $system_in does not exist"
		continue
	fi

	if [ ! -d $test_dir ]; then
		# Make a directory for the test
		mkdir $test_dir
	fi

	if [ ! -d $test_dir ]; then
		echo "Failed to create directory $test, skipping"
		continue
	fi

	# Copy over the necessary files
	cp $reference_dir/${test_dir}/${params_file} $test_dir
	cp $reference_dir/${test_dir}/${system_in} $test_dir

	# Move into that directory
	cd $test_dir

	for i in ${num_threads[@]}
	do

		echo -e "\tProcessing OMP_NUM_THREADS=$i"
		export OMP_NUM_THREADS=$i

		# Run MOLDY and append stdout and stderr to file
		../$moldy_exe &> $output_file

		# Check that the files agree to within a given number of significant digits
		system_out_diff="`$compare system.out ../$reference_dir/$test_dir/system.out`"
		# Diff output66 but ignore the CPU-TIME used as this obviously varies
		output_66_diff="`$compare -i "CPU-TIME USED" output66.txt ../$reference_dir/$test_dir/output66.txt`"

		if [ "$system_out_diff" == '' -a "$output_66_diff" == '' ]; then
			echo -e "\tTest PASSED"
		else
			echo -e "\tTest FAILED"
			fail_dir="P=${i}FAILED"
			mkdir "$fail_dir"
			cp system.out $fail_dir
			cp output66.txt $fail_dir
			cp run.out $fail_dir
			echo $system_out_diff > $fail_dir/system.out.diff
			echo $output_66_diff > $fail_dir/output66.diff
			echo -e "\tRun output saved to $test_dir/$fail_dir"
		fi

		rm system.out
		rm output66.txt

		echo -e "\tFinished OMP_NUM_THREADS=$i"
	done

	cd ..

	echo "Finished test $test_dir"
done
