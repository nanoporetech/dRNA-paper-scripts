#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from gtfparse import read_gtf_as_dataframe
import pandas as pd
import numpy as np
from collections import defaultdict
from wub.vis import report

"""
Parse command line arguments.
"""
parser = argparse.ArgumentParser(
    description='Plot normalised exon coverages.')
parser.add_argument(
    '-r', metavar='out_pdf', type=str, help="Output pdf.", default="gtf_plot_exon_cov.pdf")
parser.add_argument('gtf', metavar='gtf', type=str, help="Input gtf file.")


def plot_arrays(self, data_map, title="", xlab="", ylab="", marker='.', legend_loc='best', legend=True):
    """Plot multiple pairs of data arrays.

    :param self: object.
    :param data_map: A dictionary with labels as keys and tupples of data arrays (x,y) as values.
    :param title: Figure title.
    :param xlab: X axis label.
    :param ylab: Y axis label.
    :param marker: Marker passed to the plot function.
    :param legend_loc: Location of legend.
    :param legend: Plot legend if True
    :returns: None
    :rtype: object
    """
    fig = self.plt.figure()

    cmap = self.plt.cm.tab20
    color = iter(cmap(np.linspace(0, 1, len(data_map))))
    for label, data_arrays in data_map.items():
        self.plt.plot(data_arrays[0], data_arrays[1], marker, label=label, color=next(color))

    if legend:
        self.plt.legend(loc=legend_loc)
    self._set_properties_and_close(fig, title, xlab, ylab)


if __name__ == '__main__':
    args = parser.parse_args()

    df = read_gtf_as_dataframe(args.gtf)
    df_trs = df[df["feature"] == "exon"]

    et = defaultdict(lambda: defaultdict(dict))

    for row in df_trs.itertuples():
        et[row.gene_id][row.transcript_id][int(row.exon_number)] = float(row.cov)

    etp = defaultdict(lambda: defaultdict(list))
    for gene, trs_info in sorted(et.items(), key=lambda x: x[0]):
        for trs, exon_info in trs_info.items():
            exon_numbers = np.array(list(exon_info.keys()), dtype=int)
            exon_cov = np.array(list(exon_info.values()), dtype=float)
            # exon_cov = exon_cov / np.sum(exon_cov)
            exon_cov = np.log(exon_cov + 1)
            etp[gene][trs].extend([exon_numbers, exon_cov])

    plotter = report.Report(args.r)

    for gene, data_map in etp.items():
        plot_arrays(plotter, data_map, title=gene, xlab='Exon number', ylab='log(average base coverage+1)',
                    marker='o-', legend_loc='best', legend=True)

    plotter.close()
