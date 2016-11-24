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

def print_stats(data, n_col):
	data = data[:,n_col]
	mu = data.mean()
	sigma = data.var()
	N = data.shape[0]
	conf_bottom = np.percentile(data, 5)
	conf_top = np.percentile(data, 95)
	print "%.6f	%.6f	%.6f	%.6f	%d" % (mu, np.sqrt(sigma), conf_bottom, conf_top, N)

# Setting contants
calculate_data = 0 # If zero, calculate. Otherwise, read in.

# Reading arguments
try:
	#print sys.argv
	geno_filename = str(sys.argv[1])
	anc_filename = str(sys.argv[2])
	prefix = str(sys.argv[3])
	n_col_ulchihg = int(sys.argv[4])
	n_bootstrap = int(sys.argv[5])	
except ValueError:
	raise MyError, "wrong number of arguments, need 4"

mean_filename = "%s/prop_means.txt" % prefix
median_filename = "%s/prop_medians.txt" % prefix
props_filename = "%s/props.txt" % prefix

# Reading genotype file
with open(geno_filename) as f:
	ncols = len(f.readline().split(' '))

geno_data_raw = np.loadtxt(geno_filename, delimiter=' ', usecols=range(6,ncols))
n_snps = geno_data_raw.shape[1]/2
n_geno = geno_data_raw.shape[0]
print >> sys.stderr, "%d individuals" % n_geno
geno_data = np.empty([n_snps,n_geno])
for i_snp in range(0, n_snps):
	geno_data[i_snp,:] = np.add(geno_data_raw[:,i_snp*2]-1, geno_data_raw[:,i_snp*2+1]-1) # In {0,1,2}

# Reading ancestral population frequencies
anc_data = np.loadtxt(anc_filename, delimiter = ' ')
n_anc = anc_data.shape[1]
print >> sys.stderr, "%d ancestral components" % n_anc

# Bootstrapping
bootstrap_means = np.empty([n_bootstrap, n_anc])
bootstrap_medians = np.empty([n_bootstrap, n_anc])
bootstrap_props = np.empty([n_bootstrap, n_geno])

for i_bootstrap in range(0, n_bootstrap):
	print >> sys.stderr, "Bootstrap %d" % i_bootstrap
	
	out_filename = "%s/Props/props%d.txt" % (prefix, i_bootstrap)
	if (calculate_data == 0):
		random_indices = np.random.randint(n_snps, size=n_snps)
		random_indivs = np.random.randint(n_geno, size=n_geno)
		props_minim = estimate_proportions(geno_data[random_indices,:][:,random_indivs], anc_data[random_indices,:])
		np.savetxt(out_filename, props_minim, fmt='%.6f')
	else:
		props_minim = np.loadtxt(out_filename, delimiter=' ')
	
	bootstrap_means[i_bootstrap,:] = np.mean(props_minim, axis=0) # Mean of each column
	bootstrap_medians[i_bootstrap,:] = np.median(props_minim, axis=0) # Median of each column
	bootstrap_props[i_bootstrap,:] = props_minim[:,n_col_ulchihg]

# Save means and medians
np.savetxt(mean_filename, bootstrap_means, fmt='%.6f')
np.savetxt(median_filename, bootstrap_medians, fmt='%.6f')
np.savetxt(props_filename, bootstrap_props, fmt='%.6f')

# Compute confidence intervals
print_stats(bootstrap_means, n_col_ulchihg)
print_stats(bootstrap_medians, n_col_ulchihg)
