#!/bin/bash
# Merge several Plink files
# Need normal Plink, outputs binary plink
# First argument is output, rest is files to be merged

if [ $# -lt 3 ]
then
	>&2 echo "Not enough files to merge, need at least 2"
	>&2 echo "Note: first argument is output file"
	exit 2
fi

INTER=$(mktemp temp.XXXXXX)
touch "${INTER}.bed"
touch "${INTER}.bim"
touch "${INTER}.fam"
touch "${1}_mergescript.log"

echo "Merging ${2} and ${3}"
plink --bfile $2 --bmerge $3.bed $3.bim $3.fam --make-bed --out ${INTER} >"${1}_mergescript.log" 2>/dev/null
# If there are triallelic SNPs, exclude them and then merge
if [ -a ${INTER}.missnp ]
then
	N_TRI=$(wc -l ${INTER}.missnp | awk '{print $1}')
	echo "Dealing with ${N_TRI} triallelic sites in ${2} and ${3}"
	cat ${INTER}.missnp > ${1}_triallelic.txt
	plink --bfile $2 --exclude ${INTER}.missnp --make-bed --out ${2}_nomis >>"${1}_mergescript.log"
	plink --bfile $3 --exclude ${INTER}.missnp --make-bed --out ${3}_nomis >>"${1}_mergescript.log"
	plink --bfile ${2}_nomis --bmerge ${3}_nomis.bed ${3}_nomis.bim ${3}_nomis.fam --make-bed --out ${INTER} >>"${1}_mergescript.log"
	
	rm ${2}_nomis*
	rm ${3}_nomis*
fi
# Move intermediate file back to final result
mv "${INTER}.bed" "${1}.bed"
mv "${INTER}.bim" "${1}.bim"
mv "${INTER}.fam" "${1}.fam"

for num in ${@:4}
do
	echo "Merging ${num}"
	plink --bfile ${1} --bmerge ${num}.bed ${num}.bim ${num}.fam --make-bed --out ${INTER} >>"${1}_mergescript.log" 2>/dev/null
	# If there are triallelic SNPs, exclude them and then merge
	if [ -a ${INTER}.missnp ]
	then
		N_TRI=$(wc -l ${INTER}.missnp | awk '{print $1}')
		echo "Dealing with ${N_TRI} triallelic sites in ${num}"
		cat ${INTER}.missnp >> ${1}_triallelic.txt
		plink --bfile $1 --exclude ${INTER}.missnp --make-bed --out ${1}_nomis >>"${1}_mergescript.log"
		plink --bfile $num --exclude ${INTER}.missnp --make-bed --out ${num}_nomis >>"${1}_mergescript.log"
		plink --bfile ${1}_nomis --bmerge ${num}_nomis.bed ${num}_nomis.bim ${num}_nomis.fam --make-bed --out ${INTER} >>"${1}_mergescript.log"

		rm ${1}_nomis*
		rm ${num}_nomis*
	fi
	# Move intermediate file back to final result
	mv "${INTER}.bed" "${1}.bed"
	mv "${INTER}.bim" "${1}.bim"
	mv "${INTER}.fam" "${1}.fam"
done

rm ${INTER}*
