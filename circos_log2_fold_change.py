#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import argparse
import numpy as np
import pandas as pd

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Calculate log2 fold change from two circos files.')
parser.add_argument(
    '-a', metavar='input_a', type=str, help="Input circos file A.", required=True)
parser.add_argument(
    '-b', metavar='input_b', type=str, help="Input circos file B.", required=True)


if __name__ == '__main__':
    args = parser.parse_args()

    fields = ["chrom", "start", "end", "count"]
    pseudocount = 10**-7

    dfa = pd.read_csv(args.a, sep="\t", names=fields)
    dfb = pd.read_csv(args.b, sep="\t", names=fields)

    counts_a = np.array(dfa["count"], dtype=float)
    counts_b = np.array(dfb["count"], dtype=float)

    counts_a = counts_a / np.sum(counts_a)
    counts_b = counts_b / np.sum(counts_b)

    counts_a = counts_a + pseudocount
    counts_b = counts_b + pseudocount
    ratio = counts_a / counts_b

    ratio = np.log2(ratio)

    dfa["count"] = ratio

    dfa.to_csv("processed_ratio.tsv", sep="\t", index=False)
