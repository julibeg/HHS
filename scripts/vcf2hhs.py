#!/usr/bin/env python3
"""
Accepts the stream of a *homozygous* VCF file as input and extracs the variant
matrix, which is printed to STDOUT. The output has all sample IDs in the
first line while every consecutive line has the following format:
`CHR_POS_REF_ALT:000001000010011...`
with `0` for the reference, `1` for an alternative variant, and `2` for a
missing genotype.
Usually, this would be used in conjunction with `bcftools` (e.g.
`bcftools convert -Ov file.bcf.gz | vcf2hhs.py | gzip -c > file.hhs.gt.gz`).
For large datasets (100k+  variants and 20k+ samples) it might take up to half
an hour.
"""

import sys
eprint = sys.stderr.write


def parse_line(line):
    line = line.strip().split()
    variant = '_'.join(line[:2] + line[3:5])
    gt_line = ''
    for gt_field in line[9:]:
        gt = gt_field.split(':')[0]
        # we assume all variants are homozygous!
        if gt[0] == '0':
            gt_line += '0'
        elif gt[0] == '.':
            gt_line += '2'
        else:
            gt_line += '1'
    # write parsed line
    print(f'{variant}:{gt_line}')


with sys.stdin as input:
    # skip header
    for line in input:
        if line.startswith('#'):
            old_line = line
        else:
            break

    # the last line of the header holds the sample IDs
    samples = old_line.strip().split()[9:]
    eprint(f'header processed -- {len(samples)} samples\n')
    print(','.join(samples))

    # process the VCF entries
    parse_line(line)
    counter = 1
    for line in input:
        parse_line(line)
        counter += 1
        eprint('\033[2K\033[1G')
        eprint(f'{counter} entries processed')
        sys.stderr.flush()

eprint('\ndone\n')
