# Installation
Python 3.x required
plus the following libraries:
numpy, scipy, pandas, tables, IPython

# Usage
## 1. Permutations

### arguments: 
num_cores num_perms batch_size chunk_size input_csv output_h5_file

### example usage:
```bash
python3 permutations.py 20 10000 200 100000 450k_Mvalue_matrix.csv permutations.h5
```
### notes: 
* batch_size is number of permutations to run per core
* batch_size must be exactly divisible by num_cores (% == 0)
* num_perms must be exactly divisible by batch_size (% == 0) - recommended to leave as 10000
* chunk_size is number of rows of input data to read at a time : best performance when reading in all rows, but this requires approx (4 * size of input csv) of free system memory. If encountering memory related problems, try setting this to a smaller number

The values given in the example above have been found to perform well on a Intel i7 system with 32 cores and 16GB memory, with an input file of up to 4GB in size read in batches of (~ 0.25 * number of rows)  - i.e. 25% at a time.

## 2. Subsampling
### arguments: 
num_cores input_h5 n_cases start end reps output_csv

* start : start hd5 node - leave as 0
* end : end hd5 node - leave as 50
* reps: number of times to repeat subsampling - leave as 100

### example usage:
```bash
python2 subsampling.py 10 permutations.h5 0 50 100 results.csv
```
