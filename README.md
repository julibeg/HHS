# HHS

Hungry-Hungry SNPos (HHS) -- a "cannibalistic" elimination algorithm for finding resistance-associated genetic variants 
from a binary genotype matrix (must be homozygous, e.g. bacterial). This approach is more robust towards co-occurring 
resistance than many correlation-based methods (e.g. GWAS or machine learning techniques). It needs three input files:
the genotype matrix, the phenotypes, and a symmetric matrix of pairwise distances.

The file holding the genotype matrix should look like
```
SampleID1,SampleID2,SampleID3,...
VarID1:10X0101X...
VarID2:X00X1001...
...
```

with the sample names in the first line and the genotypes of the single variants below. Any character other than '0'
or '1' can represent missing values. 

The phenotypes should be provided in a CSV file with '0' and '1' denoting the 
absence or presence of the phenotype, respectively. Missing values are not allowed in the phenotype file and such samples
should be removed from all files before running the program. 

To account for population structure, a symmetric pairwise distance matrix of the samples should be provided 
in a CSV file with the sample names in the header (i.e. the first row) and in the index column (and '0's in the main diagonal).

## Installation

For many 64-bit Linux systems the binary in `/bin` should work right out of the box. Otherwise, please install Rust,
clone the repository and run `cargo build --release` in the root directory. The generated binary will be
in `/target/release`


## Usage example
Input files to test the installation are in `example_files.tar.gz`. They are from a highly multidrug-resistant *M. tuberculosis* dataset with phenotype data for resistance to isoniazid.
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
indices yields `1674263, 2155168, 2155169, 2155259, 2155541, 2155607, 2155786`. These genomic positions are all 
within the *inhA* and *katG* genes, which makes sense as isoniazid attacks the gene product of *inhA*
and is activated by the product of *katG*.
