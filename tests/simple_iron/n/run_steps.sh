#!/bin/bash

tools_reldir="../../../tools/"	# Path to directory where all the tools binaries are stored
moldy_exe="../moldy.iron"	# The moldy binary to use for running the sim, must be compiled for right potential
preproc_file="preprocess.sh"	# Preprocess file that can be executed to perform step specific actions
input_dir="../system.in"	# Directory with system input files

my_wd="`pwd`"

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

echo "Binaries directory is $bin_dir"

input_files=$input_dir/*

for input in $input_files
do
		echo "Processing input $input"

		if [ -a $preproc_file ]; then
			echo -e "\tPreprocessing"
			./$preproc_file
			echo -e "\tFinished preprocessing"
		fi

		# Let's make sure the number of atoms in system.in and params.in agree
		java -classpath "$bin_reldir" FixParams $input_dir/$input params.in

		# Run MOLDY
		$moldy_exe
		
		echo "Finished input $input"
i=$[$i+1]	
done
