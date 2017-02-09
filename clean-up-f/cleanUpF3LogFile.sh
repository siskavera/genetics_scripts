#!/bin/bash
# Cleans up output logfile from AdmixTools
# $1 logfile name

printf "Source1\tSource2\tTarget\tf3\tSE\tZ\tSNPs\n"
grep "result:" "${1}" | awk 'BEGIN{OFS = "\t"}{print $2,$3,$4,$5,$6,$7,$8}'
