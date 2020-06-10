# HHS

Hungry-Hungry SNPos (HHS) - a "cannibalistic" elimination algorithm for finding resistance-associated genetic variants 
from a binary genotype matrix. This approach is more robust towards co-occurring resistance than most 
correlation-based methods (e.g. GWAS or machine learning techniques). 

The file holding the genotype matrix should look like
```
10X0101X 
X00X1001
```

with one SNP per row and one sample per column. `X` denotes missing values (any character other than `0` or `1` can be
selected to represent NAs). The file can also be transposed (one sample per row and one SNP per column).

The file with the phenotypes should have the same layout but feature only a single line. Here, no missing data are 
allowed (remove samples with missing phenotypes before analysis).

To account for population structure a symmetric pairwise distance matrix of the samples should be provided 
in a `.csv` file with `0`s in the main diagonal. The distances can be obtained from a phylogeny or calculated 
from the genotype matrix with https://github.com/julibeg/binary-hamming-dist.

## Installation

For some Linux systems the binary in `/bin` might work right out of the box. Otherwise, please install Rust,
clone the repository and run `cargo build --release` in the root directory. The generated binary should be
in `/target/release`


## Usage example
Input files to test the installation are in `example_files.tar.gz`. They are from a *M. tuberculosis* dataset with resistance 
to isoniazid.
The pairwise distances are already provided in `inh.dist`, but could be generated with 
```
./binary_hamming_dist -i inh.gt -n 2 -o inh.dists -t N
```
using [binary_hamming_dist](https://github.com/julibeg/binary-hamming-dist) with `N` threads. 
For running 30,000 iterations (the default) of HHS please type 
```
./hhs -s inh.gt -p inh.phen -d inh.dists -n 2 -t N --p1g1_filter 3 -T
```
again with `N` denoting the desired number of threads used. 
`-T` is required because `inh.gt` holds a SNP per column instead of per row, which is the default.
The program should generate `16329, 20794, 20795, 20805, 20815, 20819, 20836` which are the indices of the final
SNPs remaining. The SNP positions of the dataset can be found in `inh.gt.pos`. Looking up the respective
indices yields `1674263, 2155168, 2155169, 2155259, 2155541, 2155607, 2155786`. These positions are all 
within the *inhA* and *katG* genes, which makes sense as isoniazid attacks the gene product of *inhA*
and is activated by the product of *katG*.
