# ADMIXTURE projection scripts

Projects samples on pre-defined ADMIXTURE components by numerically optimising the likelihood function from (1), using the logic of (2) and (3).

## Project single sample
### Usage:
project_admixture.py \<genotypes.ped\> \<ancestral_components.P\> \<output_prefix\>

### Inputs:
* \<genotypes.ped\>: 12-encoded PED file. Projection will be calculated for all samples in the file. Reference bases and SNP locations need to be exactly the same as in the file used in ADMIXTURE!
* \<ancestral components.P\>: Output from ADMIXTURE, containing allele frequencies for the ancestral components.
* \<output_prefix\>: Prefix for output.

### Output:
* \<prefix\>_props.txt: Contains inferred ancestry proportions for all samples (rows) and components (columns).

## Bootstrap sample(s), resampling SNPs only

### Usage:
bootstrap_snpsonly.py \<genotypes.ped\> \<ancestral_components.P\> \<output_prefix\> \<n_replicates\>

### Inputs:
* \<genotypes.ped\>: 12-encoded PED file. Projection will be calculated for all samples in the file. Reference bases and SNP locations need to be exactly the same as in the file used in ADMIXTURE!
* \<ancestral components.P\>: Output from ADMIXTURE, containing allele frequencies for the ancestral components.
* \<output_prefix\>: Prefix for output.
* \<n_replicates\>: Number of bootstrap replicates.

### Output:
* \<prefix\>_props_i.txt: Contains inferred ancestry proportions for all samples (rows) and components (columns) from replicate i.
* \<prefix\>_summary.txt: Contains summary statistics (mean, standard deviation, 95% confidence interval bottom and top) per sample and component. Tab-separated, with header.

### Requirements:
* python 2 (I used 2.7.6)
* numpy
* scipy

### References
1. D. H. Alexander, J. Novembre, K. Lange, Fast model-based estimation of ancestry in unrelated individuals. Genome Res. (2009), doi:10.1101/gr.094052.109.
2. M. E. Allentoft et al., Population genomics of Bronze Age Eurasia. Nature. 522, 167–172 (2015).
3. M. Sikora et al., Population Genomic Analysis of Ancient and Modern Genomes Yields New Insights into the Genetic Ancestry of the Tyrolean Iceman and the Genetic Structure of Europe. PLOS Genet. 10, e1004353 (2014).

### Note:
bootstrap_samples.py: Used to project a random selection of samples multiple times, on a random selection of SNP-s. Not cleaned up
