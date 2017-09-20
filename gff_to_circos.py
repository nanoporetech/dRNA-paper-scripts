#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import argparse
import numpy as np

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Convert BED file produced by bedtools coverage to circos format.')
parser.add_argument(
    '-i', metavar='input', type=str, help="Input bed.", required=True)
parser.add_argument(
    '-x', action="store_true", default=False, help="Do not take log.")



if __name__ == '__main__':
    args = parser.parse_args()

    fh = open(args.i, 'r')

    for line in fh:
        records = line.split("\t")
        chrom = "chr" + records[0]
        if chrom == "chrMito":
            chrom = "chrM"
        start = int(records[3]) - 1
        end = int(records[4]) - 1
        if not args.x:
            count = np.log(int(records[9]) + 1.0)
        else:
            count = int(records[9])


        print("{}\t{}\t{}\t{}".format(chrom, start, end, count))
