#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from gtfparse import read_gtf_as_dataframe
import pandas as pd
import numpy as np

"""
Parse command line arguments.
"""
parser = argparse.ArgumentParser(
    description='Convert stringtie GTF into tab separated counts file.')
parser.add_argument(
    '-t', metavar='out_tsv', type=str, help="Output tab separated file.", default="gtf_to_counts.tsv")
parser.add_argument('gtf', metavar='gtf', type=str, help="Input gtf file.")


if __name__ == '__main__':
    args = parser.parse_args()

    df = read_gtf_as_dataframe(args.gtf)
    df_trs = df[df["feature"] == "transcript"]

    df_trs = pd.DataFrame(df_trs, columns=["transcript_id", "FPKM"])

    df_trs["Reference"] = df_trs["transcript_id"]
    df_trs["Count"] = np.array(df_trs["FPKM"], dtype=float)
    df_trs = pd.DataFrame(df_trs, columns=["Reference", "Count"])

    df_trs.sort_values(['Count'], ascending=False, inplace=True)
    df_trs.to_csv(args.t, sep="\t", index=False)
