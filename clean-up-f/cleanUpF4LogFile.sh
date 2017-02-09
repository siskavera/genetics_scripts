#!/bin/bash
# Cleans up output logfile from AdmixTools qpF4ratio
# $1 logfile name

printf "X\tA\tB\tC\tO\talpha\tSE\tZ\n"
grep "result:" "${1}" | awk 'BEGIN{OFS = "\t"}{print $4, $2, $9, $5, $3, $11, $12, $13}'
