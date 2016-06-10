#!/bin/bash
#
# Linux shell script to convert VCF to Plink format
# Output haploid (represented as fully homozygote in Plink)
# "Haploidification" by randomly choosing an allele
# $1 Input vcf
# $2 Output file name
# $3 Extra input if need reference as well

# Note to self: triallelic loci are set to missing

# Check if DP4 field has been added
if [ $(grep -c "FORMAT=<ID=DP4" $1) -eq 0 ]
then
	echo "vcf file does not contain per sample DP4 field"
	echo "Can be added using -t DP4 flag in samtools during calling"
	exit
fi

# Temporary file
TEMPFILE=$(mktemp temp.XXXXXXXXXX.ped)

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
	if (NPOPS < 1)
	{
		print "No populations in vcf file" > "/dev/stderr"
		exit 1
	}
	for (i=1; i<NPOPS; i++)
		printf "%s", $(i+9) OFS;
	if (NF) printf "%s",$NF; 
	printf ORS

	for (i=1; i<NPOPS; i++)
		printf "%s", $(i+9) OFS;
	if (NF) printf "%s",$NF; 
	printf ORS

	for (i=0; i<3; i++)
	{
		for (j=1; j<NPOPS; j++)
			printf "%s", 0 OFS;
		printf "%s", 0; printf ORS;
	}

	for (j=1; j<NPOPS; j++)
		printf "%s", 1 OFS;
	printf "%s", 1; printf ORS;
}

NR>nlines {
	REF = $4;
	ALT = $5;
	if (length(ALT) == 1) # All reads reference or biallelic site
	{
		for (i=1; i<=NPOPS; i++)
		{
			split($(i+9), array, ":")
			GT = array[nGT]
			split(array[nDP4], DP4, ",") # DP4: Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases
			nReads = DP4[1] + DP4[2] + DP4[3] + DP4[4]
			if ( nReads == 0 )
			{
				GT_ALL[i] = "0"
			}
			else if ( GT == "0/0")
			{
				GT_ALL[i] = REF;
			} 
			else if ( GT == "0/1" )
			{
				if ( rand() < 0.5 )
					GT_ALL[i] = REF
				else
					GT_ALL[i] = ALT
			} 
			else if ( GT == "1/1" )
			{
				GT_ALL[i] = ALT;
			}
			else if ( GT == "1/0" )
			{
				if ( rand() < 0.5 )
					GT_ALL[i] = REF
				else
					GT_ALL[i] = ALT
			}
			else
			{
				GT_ALL[i] = "0";
			}
		}		
			
		for (i=1; i<NPOPS; i++)
		{
			printf "%s", GT_ALL[i] OFS
		}
		if (NPOPS) printf "%s", GT_ALL[NPOPS];
		printf ORS;
		
		for (i=1; i<NPOPS; i++)
		{
			printf "%s", GT_ALL[i] OFS
		}
		if (NPOPS) printf "%s", GT_ALL[NPOPS];
		printf ORS;
	}
	else
	{
		for (i=1; i<NPOPS; i++)
			printf "%s", 0 OFS
		printf "%s", 0; printf ORS;
		
		for (i=1; i<NPOPS; i++)
			printf "%s", 0 OFS
		printf "%s", 0; printf ORS;
	}
}' $1 > "${TEMPFILE}"

# Transpose
ncols=`head -n 1 "${TEMPFILE}" | wc -w` 
for (( i=1; i <= $ncols; i++))
do
  awk '{printf ("%s%s", tab, $'$i'); tab="\t"} END {printf ORS}' "${TEMPFILE}"
done > $2.ped
rm "${TEMPFILE}"

# Write map file
awk -v nlines=$nlines '
BEGIN{OFS = "\t"}
NR>nlines {nLength = length($2);
	factor = 10^(nLength - 5);
	rounded = sprintf("%.0f", $2 / factor) * factor;
	printf "%s\t%s\t%f\t%s\n", $1, $3, rounded /100000000, $2}' $1 | sed 's/^chr//' - > $2.map

# Write reference if three or more arguments
if [ $# -gt 2 ]
then
	awk -v nlines=$nlines '
	BEGIN{
		OFS = "\t";
		printf "%s", "Ref"OFS"Ref"OFS"0"OFS"0"OFS"0"OFS"1"OFS;
	}
	NR>nlines {
		REF = $4;
		printf "%s", REF OFS REF OFS;
	}
	END{
		printf "%s", ORS
	}' $1 >> $2.ped
fi
