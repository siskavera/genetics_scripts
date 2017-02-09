#!/bin/bash
# Cleans up output logfile from AdmixTools qpF4ratio
# $1 logfile name

printf "outgroup\tX\tasia_hg\tasia_modern\teu_alt\teu_hg\teu_modern\talpha\tSE\tZ\n"
grep "result:" "${1}" | awk 'BEGIN{OFS = "\t"}{print $2, $3, $4, $5, $8, $9, $10, $11, $12, $13}'
