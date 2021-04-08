#!usr/bin/env bash

gt_file=$1

if (file $gt_file | grep -q compressed )
then
     my_cat=zcat
else my_cat=cat
fi

 cat <($my_cat $gt_file | head -n1 | sed 's/^/,/') <(paste -d',' <($my_cat $gt_file | head -n1 | tr ',' '\n') <(binary_hamming_dist -i <($my_cat $gt_file | sed '1d' | cut -d: -f2) -T -t 6))