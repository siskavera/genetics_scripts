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
plink --file $2 --merge $3.ped $3.map --make-bed --out ${INTER} >"${1}_mergescript.log" 2>/dev/null
# If there are triallelic SNPs, exclude them and then merge
if [ -a ${INTER}.missnp ]
then
	echo "Dealing with triallelic sites in ${2} and ${3}"	
	cat ${INTER}.missnp > ${1}_triallelic.txt
	plink --file $2 --exclude ${INTER}.missnp --recode --out ${2}_nomis >"${1}_mergescript.log"
	plink --file $3 --exclude ${INTER}.missnp --recode --out ${3}_nomis >"${1}_mergescript.log"
	plink --file ${2}_nomis --merge ${3}_nomis.ped ${3}_nomis.map --make-bed --out ${INTER} >"${1}_mergescript.log"

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
	plink --bfile ${1} --merge ${num}.ped ${num}.map --make-bed --out ${INTER} >"${1}_mergescript.log" 2>/dev/null
	# If there are triallelic SNPs, exclude them and then merge
	if [ -a ${INTER}.missnp ]
	then
		echo "Dealing with triallelic sites in ${num}"
		cat ${INTER}.missnp > ${1}_triallelic.txt
		plink --file $1 --exclude ${INTER}.missnp --recode --out ${1}_nomis >"${1}_mergescript.log"
		plink --file $num --exclude ${INTER}.missnp --recode --out ${num}_nomis >"${1}_mergescript.log"
		plink --file ${1}_nomis --merge ${num}_nomis.ped ${num}_nomis.map --make-bed --out ${INTER} >"${1}_mergescript.log"

		rm ${1}_nomis*
		rm ${num}_nomis*
	fi
	# Move intermediate file back to final result
	mv "${INTER}.bed" "${1}.bed"
	mv "${INTER}.bim" "${1}.bim"
	mv "${INTER}.fam" "${1}.fam"
done

rm ${INTER}*
