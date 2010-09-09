#!/bin/bash

if [ ! $# == 3 ]; then
	echo "Usage: $0 params_file parameter_token new_value"
	exit
fi

# Get the command line parameters
params_file=$1
param_token=$2
param_value=$3

if [ ! -f $params_file ]; then
	echo "Params file '$params_file' does not exist"
	exit
fi

old_value=`grep -m 1 "$param_token" $params_file | grep -Eo "=[[:blank:]]*[^[:blank:]]+"`

echo "Changing $param_token $old_value to $param_value"

if [ -n "$old_value" ]; then
	sed -e "s/\(nbrupdate[[:blank:]]*=[[:blank:]]*\)[^[:blank:]]*\([[:blank:]]*#*.*\)/\1${param_value}\2/" $params_file > tmpfile && mv tmpfile $params_file
else
	# The token wasn't found so just add the parameter at the end of the file
	echo -e "$param_token\t= $param_value" >> $params_file
fi
