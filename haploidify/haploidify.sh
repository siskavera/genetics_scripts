#!/bin/bash
#
# Linux shell script to haploidify ped file
# Plink treates completely homozygous files as haploid
# Note: need to remove carriage return from windows files

awk -v seed=$RANDOM 'BEGIN{
	OFS="\t";
	ORS="\n";
	srand(seed);
}
{
# Deleting carriage return in case non-Unix system was used
	gsub(/\r/,"",$NF)
# Printing sample info
	for (i=1; i<=6; i++)
		printf "%s", $i OFS
# Printing haploidified data
	for (i=7; i<=(NF-2); i += 2)
	{
		if ( rand() < 0.5 )
			printf "%s", $i OFS $i OFS
		else
			printf "%s", $(i+1) OFS $(i+1) OFS
	}
# Last record
	if ( rand() < 0.5 )
		printf "%s", $(NF-1) OFS $(NF-1) ORS
	else
		printf "%s", $(NF) OFS $(NF) ORS
}' $1