#!/bin/bash
# Bash script to deaminate a .ped file
# Change T -> C
# Change G -> A
# Delimiter !!!!!!!!!!

# Changes in sequence file
# Need spaces before and after so that sample IDs don't get matched
# Need to run the command twice to match both alleles (a match eats up the space)
# End of line!
sed -e 's/ T\( \|$\)/ C\1/g' -e 's/ G\( \|$\)/ A\1/g' -e 's/ T\( \|$\)/ C\1/g' -e 's/ G\( \|$\)/ A\1/g' "${1}"

