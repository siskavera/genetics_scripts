#!/bin/bash
#
# Linux shell script to count when random() matters

# Note to self: triallelic loci are set to missing

# Check if DP4 field has been added
if [ $(grep -c "FORMAT=<ID=DP4" $1) -eq 0 ]
then
	echo "vcf file does not contain per sample DP4 field"
	echo "Can be added using -t DP4 flag in samtools during calling"
	exit
fi

# Finding last line of header
nlines=$(grep -n "^#" $1 | tail -1 | sed 's/:.*//')

# Finding #field for GT and DP4 in INFO
nGT=$(grep -v -m 1 "^#" $1 | awk '{print $9}' | sed 's/:/\n/g' | grep -n "GT" | cut -d ":" -f 1)
nDP4=$(grep -v -m 1 "^#" $1 | awk '{print $9}' | sed 's/:/\n/g' | grep -n "DP4" | cut -d ":" -f 1)

input="${nlines} ${nGT} ${nDP4}"

# Conversion to transposed Ped
awk -v input="${input}" '
BEGIN {
	split(input, array, " ");
	nlines=array[1];
	nGT=array[2];
	nDP4=array[3];
}	
NR==nlines {
	NPOPS = NF - 9
}
NR>nlines {
	REF = $4;
	ALT = $5;
	if (length(ALT) == 1) # All reads reference or biallelic site
	{
		for (i=1; i<=NPOPS; i++)
		{
			# Get read counts
			split($(i+9), array, ":")
			split(array[nDP4], DP4, ",") # DP4: Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases
			nREF = DP4[1] + DP4[2]
			nALT = DP4[3] + DP4[4]
			pALT = nALT / (nREF + nALT)

			if ( (nREF + nALT) == 0 ) # If no reads, set to missing
			{
			}
			else if ( pALT > 0 && pALT < 1 )
			{
				printf "%f\n", pALT
			}
		}		

	}

}' $1

