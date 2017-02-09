#!/usr/bin/python

import numpy as np
from scipy.optimize import minimize
from scipy import stats
import sys

def fun(x, f, g):
	x = x/np.sum(x)
	return -(g.dot(np.log(f.dot(x)))+(2-g).dot(np.log((1-f).dot(x))))

def estimate_proportions(geno_data, anc_data, use_random_indices = False):
	# Get contants
	n_geno = geno_data.shape[1]
	n_anc = anc_data.shape[1]

	# Random sampling
	if use_random_indices:
		n_snps = geno_data.shape[0]
		random_indices = np.random.randint(n_snps, size=n_snps)
		geno_data = geno_data[random_indices,:]
		anc_data = anc_data[random_indices,:]
	
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
	n_replicates = int(sys.argv[4])
except IndexError:
	sys.exit("Usage: project_admixture.py <genotypes.ped> <ancestral_components.P> <output_prefix> <n_replicates>")


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
all_props = np.empty([n_geno, n_anc, n_replicates])
print >> sys.stderr, "Starting projections, %d replicates" % n_replicates
for i_rep in range(n_replicates):
	print >> sys.stderr, "%d/%d" % (i_rep+1, n_replicates)
	out_filename = "%s_props_%d.txt" % (prefix, i_rep)
	props_minim = estimate_proportions(geno_data, anc_data, True)
	all_props[:,:,i_rep] = props_minim
	np.savetxt(out_filename, props_minim, fmt='%.6f')

print >> sys.stderr, "calculating summary statistics"
summary_filename = "%s_summary.txt" % (prefix)
with open(summary_filename, "w") as outfile:
	print >> outfile, "ind\tcomponent\tmean\tstd\t95bottom\t95top"
	for i_ind in range(n_geno):
		for i_comp in range(n_anc):
			this_data = all_props[i_ind,i_comp,:].reshape((n_replicates,1))
			this_mean = np.mean(this_data)
			this_std = np.std(this_data)
			this_bottom = np.percentile(this_data, 2.5)
			this_top = np.percentile(this_data, 97.5)
			print >> outfile, "%d\t%d\t%.6f\t%.6f\t%.6f\t%.6f" % (i_ind, i_comp, this_mean, this_std, this_bottom, this_top)