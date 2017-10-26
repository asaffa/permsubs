# Description
Python command-line tools for performing simulation extrapolation using DNA methylation array data (M-values) 
for the purpose of estimating a genome-wide significance threshold for EWAS: http://onlinelibrary.wiley.com/doi/10.1002/gepi.22086/abstract 

1. Permutations - permute case control labels and perform t-tests across all measured probes, return absolute t values, store in HDF5 compressed file
2. Subsampling - for each permutation, sample from distribution of calculated p values for range of sampling proportions (0.1 to 1) and record the minimum p for each. Next, calculate 5% point for each density across all permutations, repeat in order to obtain mean 5% value for each density across all permutations.

# Installation
Python 3.x required
plus the following libraries:
numpy, scipy, pandas, tables, IPython

# Usage
## 1. Permutations

### arguments: 
num_cores num_perms batch_size chunk_size input_csv output_h5_file

* num_cores = processor cores to use
* num_perms = number of permutations to run
* batch_size = permutations to run per core
* chunk_size = rows of input data to read at a time
* input_csv = matrix of methylation M-values (rows = probes, cols = samples) (.csv file)
* output_h5 = HDF5 file to save results to 

### example usage:
```bash
python3 permutations.py 20 10000 200 100000 450k_Mvalue_matrix.csv permutations.h5
```
### notes: 

Tested using the values given in the example above, found to perform well on an Intel i7 system with 32 cores and 16GB memory, using an input file of up to 4GB in size read in 4 batches (0.25 * num rows at a time)

If trying different values please note:
* batch_size must be exactly divisible by num_cores (% == 0)
* num_perms must be exactly divisible by batch_size (% == 0) 
* chunk_size  : best performance when reading in all rows, but this requires approx (4 * size of input csv) of free system memory. If encountering memory related problems, try setting this to a smaller number

## 2. Subsampling
### arguments: 
num_cores input_h5 n_cases start end reps output_csv

* num_cores = processor cores to use
* input_h5 = HDF5 file containing the results from the permutation tests
* n_cases = number of subjects who are cases in the sample
* start : start hd5 node - leave as 0
* end : end hd5 node - leave as 50
* reps: number of times to repeat subsampling - leave as 100

### example usage:
```bash
python2 subsampling.py 10 permutations.h5 0 50 100 results.csv
```
