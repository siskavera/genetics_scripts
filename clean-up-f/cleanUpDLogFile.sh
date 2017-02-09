#!/bin/bash
# Cleans up output logfile from AdmixTools qpF4ratio
# $1 logfile name

printf "Pop1\tPop2\tPop3\tPop4\tD\tZ\tSNPs1\tSNPs2\tSNPs\n"
grep "result:" "${1}" | awk 'BEGIN{OFS = "\t"}{print $2, $3, $4, $5, $6, $7, $8, $9, $10}'
