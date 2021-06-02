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

with the sample names in the first line and the genotypes of the variants below. Any character other than '0'
or '1' can represent missing values. 

The phenotypes should be provided in a CSV file with '0' and '1' denoting the 
absence or presence of the phenotype, respectively. Missing values are not allowed in the phenotype file and such samples
should be removed from all files before running the program. 

To account for population structure, a symmetric pairwise distance matrix of the samples should be provided 
in a CSV file with the sample names in the header (i.e. the first row) and in the index column (and '0's in the main diagonal).

## Installation

For many 64-bit Linux systems the binary in `/bin` should work right out of the box. Otherwise, please install Rust,
clone the repository and run `cargo build --release` in the root directory. The generated binary will be
in `/target/release`.


## Usage example
Input files to test the installation are in `example_files.tar.gz`. For a quick tutorial on how these files can be generated from a BCF/VCF, see below. The SNPs are from a highly multidrug-resistant *M. tuberculosis* dataset with phenotype data for resistance to isoniazid. The data were used in the [initial publication](https://doi.org/10.1371/journal.pcbi.1008518) of HHS.

`example_files.tar.gz` contains three items: `inh.hhs.gt.gz`, `inh.dists.csv.gz`, and `inh.phen.csv`. To run the program, download the binary (or build with cargo), navigate into the directory where the archive was extracted and type 
```
./hhs -g inh.hhs.gt.gz -p inh.phen.csv -d inh.dists.csv.gz -t 6 --p1g1_filter 3 -o inh.hhs.result
```
replacing `N` with the desired number of threads. 

As we can see from the messages printed to STDOUT, the scores converged after about 5,000 iterations. 
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

## Generate input files from VCF/BCF
The genotype matrix can be obtained from a VCF/BCF with the script in `/scripts/vcf2hhs.py`. It takes the stream of an uncompressed VCF as input. 
To use BCFs, you need [bcftools](http://samtools.github.io/bcftools/bcftools.html) and run something like 
```
$ bcftools convert -Ov my_file.bcf.gz | ./vcf2hhs.py | gzip -c > my_file.hhs.gt.gz
```
The VCF that has been used to generate the files in `example_files.tar.gz` can be downloaded with 
```
$ wget https://myfiles.lshtm.ac.uk/rest/files/public/8a8c80b5771ea6840178b1e0ca853d3e -O inh.vcf.gz
```
In order to process it, type the following (we can't use `bcftools` for this because the VCF has no header except for the sample IDs)
```
$ gunzip -c inh.vcf.gz | ./vcf2hhs.py | gzip -c > inh.hhs.gt.gz
```
which might take a couple seconds. Now that we got the binary genotype matrix, let's sneak a peek to confirm the format
```
$ zcat inh.hhs.gt.gz | head -n5 | cut -c-100
ERR2512419,ERR2512432,ERR2512444,SRR6046861,ERR2512448,ERR2512616,ERR2512464,ERR2512467,SRR6046574,E
Chromosome_31_A_G:0000000000000000000000000000000000000000000000000000000000000000000000000000000000
Chromosome_64_G_C:0000000000000000000000000000000000000000000000000000000000000000000000000000000000
Chromosome_238_G_T:000000000000000000000000000000000000000000000000000000000000000000000000000000000
Chromosome_365_C_T:000000000000000000000000000000000000000000000000000000000000000000000000000000000
```
This looks like what we expect. To run HHS, we also need a pairwise distance matrix. Ideally, the distances would be obtained
from a phylogeny, but if only the genotype matrix is available, we can calculate the pairwise binary hamming distances with [binary-hamming-dist](https://github.com/julibeg/binary-hamming-dist). It expects a header-less file holding only strings of '0's and '1's as input (while any other character denotes missing values) and writes the distances as CSV without any labels. The script in `/scripts/pw_bin_hamm_dists_from_gt.sh` does the relevant data wrangling in order to add the header and index column. We call it with 
```
./pw_bin_hamm_dists_from_gt.sh inh.hhs.gt.gz | gzip -c > inh.dists.csv.gz
```
Again, let's inspect the result to see if everything went alright
```
$ zcat inh.dists.csv.gz | head -n5 | cut -c-100
,ERR2512419,ERR2512432,ERR2512444,SRR6046861,ERR2512448,ERR2512616,ERR2512464,ERR2512467,SRR6046574,
ERR2512419,0,21,27,29,28,31,36,4,6,35,27,24,9,37,61,65,24,19,28,2,25,26,10,27,30,31,27,21,27,27,22,2
ERR2512432,21,0,17,26,18,20,26,18,21,25,18,15,22,26,49,54,24,17,18,23,18,17,12,22,15,30,28,18,28,14,
ERR2512444,27,17,0,33,23,13,31,25,27,32,10,21,29,19,57,61,31,24,23,29,24,27,19,28,26,37,34,25,35,24,
SRR6046861,29,26,33,0,34,34,42,26,28,41,34,31,30,41,66,70,32,17,34,31,34,37,20,38,35,6,28,17,34,34,2
```
This looks like a well-behaved CSV with index column and header. Now we can run HHS as described above.
