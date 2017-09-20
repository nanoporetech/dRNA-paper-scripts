#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
from collections import defaultdict, OrderedDict
from wub.vis import report
import seaborn as sns
from scipy.stats import spearmanr
import numpy as np

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description='Plot correlation against mix concentrations.')
parser.add_argument(
    '-r', metavar='report_pdf', type=str, help="Output PDF.", required=True)
parser.add_argument(
    '-t', metavar='tags', type=str, help="Data tags: tag1,tag2.", required=True)
parser.add_argument(
    'counts_one', metavar='counts_one', type=str, help="Input tab separated file.")
parser.add_argument(
    'counts_two', metavar='counts_two', type=str, help="Input tab separated file.")


def get_common_keys(data_one, data_two):
    """ Get common keys. """
    set_one = set(data_one['Reference'])
    set_two = set(data_two['Reference'])
    set_one.update(set_two)
    common_keys = sorted(list(set_one))
    return common_keys


def get_count(df, tr):
    """Get counts for a transcript."""
    tmp = df.query('Reference == "{}"'.format(tr))
    if len(tmp) == 0:
        return 0.0
    return float(tmp['Count'])


def merge_data_frames(data_one, data_two, common_keys, tags):
    """Merge two data frames."""
    tags_dict = {tags[0]: data_one, tags[1]: data_two}
    res = OrderedDict([('Reference', common_keys), ('Gene', []), (tags[0], []), (tags[1], [])])
    for tr in common_keys:
        res['Gene'].append(''.join(list(tr)[:5]))
        res[tags[0]].append(get_count(data_one, tr))
        res[tags[1]].append(get_count(data_two, tr))
    return pd.DataFrame(res)


if __name__ == '__main__':
    args = parser.parse_args()

    data_one = pd.read_csv(args.counts_one, sep="\t")
    data_two = pd.read_csv(args.counts_two, sep="\t")
    tags = args.t.split(",")

    common_keys = get_common_keys(data_one, data_two)
    data_merged = merge_data_frames(data_one, data_two, common_keys, tags)

    plotter = report.Report(args.r)

    g = sns.jointplot(tags[0], tags[1], data=data_merged, stat_func=spearmanr, kind="reg")
    plotter.pages.savefig()

    data_genes = data_merged.groupby('Gene').agg(['sum'])
    data_genes = pd.DataFrame(
        OrderedDict([(tags[0], data_genes[(tags[0], 'sum')]), (tags[1], data_genes[(tags[1], 'sum')])]))

    gg = sns.jointplot(tags[0], tags[1], data=data_genes, stat_func=spearmanr, kind="reg")
    plotter.pages.savefig()

    plotter.close()
