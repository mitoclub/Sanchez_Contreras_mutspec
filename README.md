# Reproducible analysis for Sanchez-Contreras <i>et al</i>., 2022

Repository contains the necessary data and python scripts to generate nearly all of the figures and reported p-values. Reproducibility is done using Snakemake

### Requirements
Snakemake
Anaconda/miniconda

### Data Extraction
- ['compile_data.py'](compile_data.py): Reads the individual data and/or summary files and compiles them into '.csv' files used for plotting (stored in the 'data/imported_data' subdirectory).

### Statistical Analysis
- ['compute_stats.py']: Computes p-values used to establish significance and then reports them as '.csv' files for each figure or subfigure. Specific figure is indicated in the file name (stored in the 'data/stats/' subdirectory)
- ['fold_change.R']: Computes the fold change and p-values using the 'mratios' R package used in 'Figure 3'.
- ['Dunnett_test.R']: Computes the adjusted p-value using the 'DescTools' R package and used in 'Figure 6'.

### Figure Scripts
The scripts generate each of the figures in the paper. They are contained in '.py' files with names corresponding to the relevant figure.

### Miscellaneous
- ['GlobalVars_.py']: Contains global variables used across multiple figure generation and analysis scripts.
- ['HelperFuncs_.py']: Functions used for formatting data
- ['run.py']: A standalone script that will generate the figures and statistical files without Snakemake or conda.


## Howto
- Install 'Conda' or 'Miniconda'
- Install 'snakemake' and 'mamba' (['conda install snakemake mamba'])
- Clone the repository 
- Setup environment ['snakemake --cores 1 --use-conda --conda-frontend mamba --conda-prefix .snakemake -- initializeEnvs ']
- Perform reproducibile analysis (['snakemake -s snakefile --use-conda --keep-going -j 1 --conda-prefix .snakemake'])




