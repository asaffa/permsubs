# Installation
Python 2.7 required
plus the following libraries:
numpy, scipy, pandas, tables, IPython

# Usage
## 1. Permutations
###arguments: 
num_cores batch_size input_csv group_to_analyze paired_data 
equal_variance num_permutations output_h5_file

###example usage:
```bash
python permutations.py 8 200 450k_Mvalue_matrix.csv all 0 1 2000 permutations.h5
```
###notes: 
* batch_size must be exactly divisible by numcores (% == 0)
* num_permutations must be exactly divisible by batch_size (% == 0)
* paired_data not implemented - set to 0
* equal variance not implemented - set to 1

## 2. Subsampling
###arguments: 
num_cores input_h5 n_cases output_csv

###example usage:
```bash
python subsampling.py 8 permutations.h5 120 results.csv
```
