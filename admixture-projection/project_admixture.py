#!/usr/bin/python

import numpy as np
from scipy.optimize import minimize
from scipy import stats
import sys

def fun(x, f, g):
	x = x/np.sum(x)
	return -(g.dot(np.log(f.dot(x)))+(2-g).dot(np.log((1-f).dot(x))))

def estimate_proportions(geno_data, anc_data):
	# Get contants
	n_geno = geno_data.shape[1]
	n_anc = anc_data.shape[1]
	
	# Setting up for optimisation
	props_minim = np.empty([n_geno, n_anc])
	x_init = np.ones(n_anc) / n_anc
	bnds = ((1e-5, 1-1e-5),) * n_anc
	cons = ({'type': 'eq', 'fun': lambda x: np.sum(x)-1})
	
	for i_geno in range(0, n_geno):
		is_valid = geno_data[:,i_geno] >= 0
		n_nonmissing = geno_data[is_valid,i_geno].shape[0]
	
		g = geno_data[is_valid,i_geno]
		f = anc_data[is_valid,:]
		
		res = minimize(fun, x_init, args=(f,g), bounds=bnds, method="L-BFGS-B")
		props_minim[i_geno,:] = res.x / sum(res.x)
	return props_minim

# Reading arguments
try:
	#print sys.argv
	geno_filename = str(sys.argv[1])
	anc_filename = str(sys.argv[2])
	prefix = str(sys.argv[3])
except IndexError:
	sys.exit("Usage: project_admixture.py <genotypes.ped> <ancestral_components.P> <output_prefix>")


# Reading genotype file
with open(geno_filename) as f:
	ncols = len(f.readline().split(' '))

geno_data_raw = np.loadtxt(geno_filename, delimiter=' ', usecols=range(6,ncols))
n_snps = geno_data_raw.shape[1]/2
n_geno = geno_data_raw.shape[0]
print >> sys.stderr, "Projecting %d individuals" % n_geno
geno_data = np.empty([n_snps,n_geno])
for i_snp in range(0, n_snps):
	geno_data[i_snp,:] = np.add(geno_data_raw[:,i_snp*2]-1, geno_data_raw[:,i_snp*2+1]-1) # In {0,1,2}

# Reading ancestral population frequencies
anc_data = np.loadtxt(anc_filename, delimiter = ' ')
n_anc = anc_data.shape[1]
print >> sys.stderr, "on %d ancestral components" % n_anc

# Projecting all samples
out_filename = "%s_props.txt" % (prefix)
props_minim = estimate_proportions(geno_data, anc_data)
np.savetxt(out_filename, props_minim, fmt='%.6f')
print >> sys.stderr, "projection finished, results are in %s_props.txt" % (prefix)
