#!/bin/bash
# Deaminate fasta file and write to stdout

# IUPACNucleotideCode	Base	DeaminatedBase
# A	Adenine	A
# C	Cytosine	C
# G	Guanine	A
# T (or U)	Thymine (or Uracil)	C
# R	A or G	A/A=A
# Y	C or T	C/C=C
# S	G or C	A/C=M
# W	A or T	A/C=M
# K	G or T	A/C=M
# M	A or C	A/C=M
# B	C or G or T	C/A/C=M
# D	A or G or T	A/A/C=M
# H	A or C or T	A/C/C=M
# V	A or C or G	A/C/A=M
# N	any base	N
# . or -	gap	N

# Change bases
sed 's/[GR]/A/g' $1 | sed 's/[TY]/C/g' | sed 's/[SWKBDHV]/M/g'
