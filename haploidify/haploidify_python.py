#!/usr/bin/env python

import csv
import random
import sys

if (len(sys.argv) != 4):
	print "Incorrect number of arguments, need 2:"
	print "\tInput file name"
	print "\tField number"
	print "\tDelimiter"
	sys.exit(1)

if (csv.field_size_limit() < int(sys.argv[2]) ):
	csv.field_size_limit( int(sys.argv[2])*3 )

with open('new.ped', 'wb') as csvout:
	with open(str(sys.argv[1]),'rb') as tsvin:
		# need decoding for tabs
		tsvin = csv.reader(tsvin, delimiter=str(sys.argv[3]).decode("string_escape"), quoting=csv.QUOTE_NONE)
		csvout = csv.writer(csvout, delimiter='\t', quoting=csv.QUOTE_NONE)
	
		for row in tsvin:
			# Write sample info
			for i in range(0,6):
				print ('%s\t') % (row[i]) ,
			# Write data
			for i in range(6, len(row), 2):
				if ( random.random() < 0.5 ):
					print ('%s\t%s\t') % (row[i], row[i]) ,
				else:
					print ('%s\t%s\t') % (row[i+1], row[i+1]) ,
			print('\n') ,