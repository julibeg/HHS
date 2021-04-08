#!usr/bin/env bash

gt_file=$1

# figure out which cat to use based on whether the genotype 
# matrix is gzipped.
if file $gt_file | grep -q gzip
then
     my_cat=zcat
else my_cat=cat
fi

# heavy use of process substitution:
# *) we need to get the header from the gt file and add a comma to
#    the beginning for the absent label of the index column. 
#    --> <($my_cat $gt_file | head -n1 | sed 's/^/,/')
# *) `binary_hamming_dist` expects bare binary strings of '0's and '1's
#    --> we need to remove the header and the variant IDs from the gt matrix
#    --> <($my_cat $gt_file | sed '1d' | cut -d: -f2)
# *) add the samples as index column to the calculated distances with `paste`
#    --> <(paste -d',' <($my_cat $gt_file | head -n1 | tr ',' '\n') <(...)
 cat <($my_cat $gt_file | head -n1 | sed 's/^/,/') \
     <(paste -d',' <($my_cat $gt_file | head -n1 | tr ',' '\n') \
     <(binary_hamming_dist -i <($my_cat $gt_file | sed '1d' | cut -d: -f2) -T -t 6))