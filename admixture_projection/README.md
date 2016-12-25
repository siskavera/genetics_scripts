# Admixture projection script and example files

## Contains
* MOS_Korean_Panel_pruned_eacas_mos45.8.P: Output from ADMIXTURE
* project_admixture.py: script doing the projection
* refs.txt: Reference bases corresponding to the P-file
* ulchi.log: Logfile from plink 1.07, for selecting Ulchi samples
* ulchi.sh: Wrapper script to project all Ulchi samples
* ulchi.txt: List of Ulchi samples
* ulchi_fixed.log: Logfile from plink 1.07, for using the correct reference bases and recoding to 1/2
* ulchi_fixed.ped: Resulting ped file, used by projection script
* ulchi_props.txt: Output from projection script
* bootstrap_samples.py: Used to project a random selection of samples multiple times, on a random selection of SNP-s.

## Projection script
Projects samples on pre-defined ADMIXTURE components by numerically optimising the likelihood function from (1), using the logic of (2) and (3).
### Usage:
project_admixture.py <genotypes.ped> <ancestral_components.P> <output_prefix>

### Requirements:
* python 2 (I used 2.7.6)
* numpy
* scipy

### References
1. D. H. Alexander, J. Novembre, K. Lange, Fast model-based estimation of ancestry in unrelated individuals. Genome Res. (2009), doi:10.1101/gr.094052.109.
2. M. E. Allentoft et al., Population genomics of Bronze Age Eurasia. Nature. 522, 167–172 (2015).
3. M. Sikora et al., Population Genomic Analysis of Ancient and Modern Genomes Yields New Insights into the Genetic Ancestry of the Tyrolean Iceman and the Genetic Structure of Europe. PLOS Genet. 10, e1004353 (2014).
