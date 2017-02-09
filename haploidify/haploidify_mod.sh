#!/bin/bash
#
# Linux shell script to haploidify ped file
# Plink treates completely homozygous files as haploid
# Need haploidify_python.py
# Input: plink file to be haploidified

# Find out if delimiter is space or tab
head -1 "$1" | grep -q ' '
if [ $? -eq 0 ]
then
	DELIM=" "
else
	DELIM="\t"
fi

# Count number of fields
# Used to increase limit in python if needed
FIELD_NO=$(head -1 "$1" | wc -w)

# echo "Opening " "$1" " with " "${FIELD_NO}" " fields and " "${DELIM}" " as delimiter"
# Call Python script to output haploidified version to stdout
./haploidify_python.py "$1" "${FIELD_NO}" "${DELIM}"
