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
Input files to test the installation are in `example_files.tar.gz`. They are from a highly multidrug-resistant *M. tuberculosis* dataset with phenotype data for resistance to isoniazid. The data were used in the initial publication of the algorithm (https://doi.org/10.1371/journal.pcbi.1008518).

The archive `example_files.tar.gz` contains three files: `inh.hhs.gt.gz`, `inh.dists.csv.gz`, and `inh.phen.csv`.

To run 30,000 iterations (the default) of HHS, download the binary (or build with cargo), navigate into the 
directory where the archive was extracted and type 
```
./hhs -g inh.hhs.gt.gz -p inh.phen.csv -d inh.dists.csv.gz -t 6 --p1g1_filter 3 -o inh.hhs.result
```
replacing `N` with the desired number of threads. 

As we can see from the messages printed to STDOUT, the scores converged after about 10,000 iterations. 
The generated output file `inh.hhs.result` should look like 
```
./hhs -g inh.hhs.gt.gz -p inh.phen.csv -d inh.dists.csv.gz -t 6 --p1g1_filter 3 -o inh.hhs.result
VarID,score,p1g1,p0g1,dist
Chromosome_1674263_T_C,764979.25,3,0,8.498
Chromosome_2154857_C_A,14533.56,3,1,11.477
Chromosome_2155168_C_G,T,A,1311184.25,1676,12,9.818
Chromosome_2155169_T_C,1693975.75,8,0,8.798
Chromosome_2155259_C_T,A,1413699.38,3,0,10.776
Chromosome_2155541_A_G,C,1900279.00,10,0,10.793
Chromosome_2155607_C_T,3900517.75,4,0,14.937
```
with the command that has been used to invoke the program in the first line and a CSV with the results below. 
The columns in the CSV are the variant ID, final score, number of samples with the phenotype and genotype, 
number of samples with the genotype but without the phenotype, and the average pairwise distance among all samples
with the genotype. 

The variants with locations starting with `2155...` are within the *katG* gene. The other one is from *inhA*. Both genes have
been shown to be related to isoniazid resistance (the prodrug isoniazid is activated by the gene product of *katG* and 
then targets the enzyme encoded by *inhA*). Due to extensive co-occurring resistance with rifampicin, ethambutol, and streptomycin 
(among others) in this dataset, many GWAS implementations would spuriously also ascribe strong significance to variants in genes associated 
with resistance agains these drugs (like *rpoB*, *embB*, or *rrs*, respectively).
