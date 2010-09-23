#!/bin/bash

## Replace token version 1.1 						##
## by Martin Uhrin								##
## Can be called to change a setting in a params.in 			##
## Usage: ./replace_token.sh params_file parameter_token new_value 	##

if [ ! $# == 3 ]; then
	echo "Usage: $0 params_file parameter_token new_value"
	exit
fi

# Get the command line parameters
params_file=$1
param_token=$2
param_value=$3

if [ ! -f $params_file ]; then
	echo "Params file $params_file does not exist"
	exit
fi


old_value=`grep -e "$param_token[[:blank:]]*=[[:blank:]]*" $params_file`

if [ -n "$old_value" ]; then
	# Get any comments in the line
	comment=`grep -e "$param_token[[:blank:]]*=[[:blank:]]*" $params_file | grep -oe "[[:blank:]]*#.*"`

	# Now replace the value and add the comment back
	sed -e "s/\($param_token[[:blank:]]*=[[:blank:]]*\).*/\1${param_value}$comment/" $params_file > tmpfile && mv tmpfile $params_file
else
	# The token wasn't found so just add the parameter at the end of the file
	echo -e "$param_token\t= $param_value" >> $params_file
fi
