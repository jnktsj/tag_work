#!/usr/bin/env python

import os
from itertools import chain
from collections import OrderedDict, Counter, defaultdict
from argparse import ArgumentParser
import os.path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import NullFormatter

from typing import Dict
from intervaltree_bio import GenomeIntervalTree, IntervalTree, Interval


class AnnotatedInterval:
    """This class represents a half-open [,) genomic interval along with optional annotations.
    Note:
        Equality test and hashing is based on get_key() which excludes all annotations.
    """
    def __init__(self, contig, start, end, annotations = None):
        assert end >= start, "Interval end point must be greater or equal to its start point"
        self.contig = str(contig)
        self.start = int(start)
        self.end = int(end)
        self.annotations = annotations
        self._hash = hash(self.get_key())

    def get_key(self):
        return self.contig, self.start, self.end
    
    def overlap_in_bp(self, other):
        return (self.contig == other.contig) * (min(self.end, other.end) - max(self.start, other.start))

    def overlap_coord(self, other):
        return [self.contig, max(self.start, other.start), min(self.end, other.end)]

    def overlaps_with(self, other):
        return overlap_in_bp(other) > 0

    def get_midpoint(self):
        return 0.5 * (self.start + self.end)
    
    def length(self):
        return self.end - self.start
    
    def __eq__(self, other):
        return self.get_key() == other.get_key()

    def __ne__(self, other):
        return self.get_key() != other.get_key()

    def __hash__(self):
        return self._hash

    def __str__(self):
        return str(self.get_key()) if self.annotations is None else str(self.get_key()) + ': ' + str(self.annotations)

    def __repr__(self):
        return self.__str__()


class AnnotatedIntervalCollection:
    """Allows quick lookups by genomic overlap from a list of non-overlapping, annotated intervals."""
    def __init__(self, 
                 annotated_intervals_path,
                 contig_header = 'CONTIG',
                 start_header = 'START',
                 end_header = 'END',
                 is_interval_list = False):
        assert (contig_header == 'CONTIG') & (start_header == 'START') & (end_header == 'END') if is_interval_list else True
        if is_interval_list:
            annotated_intervals_df = pd.read_csv(annotated_intervals_path, sep='\t', comment='@', header=None, dtype={0: str})
            annotated_intervals_df.columns = [contig_header, start_header, end_header, 'STRAND', 'TARGET_NAME']
        else:
            annotated_intervals_df = pd.read_csv(annotated_intervals_path, sep='\t', comment='@', dtype={contig_header: str})
        self.annotated_intervals_df = annotated_intervals_df
        self.genome_interval_tree = GenomeIntervalTree()
        self.contig_dict = OrderedDict()
        for index, row in annotated_intervals_df.iterrows():
            contig = row[contig_header]
            start = row[start_header]
            end = row[end_header] + 1
            self.contig_dict[contig] = None
            if is_interval_list:
                self.genome_interval_tree.addi(contig, start, end,
                                               data=AnnotatedInterval(contig, start, end))
            else:
                annotations = row.drop([contig_header, start_header, end_header]).to_dict()
                self.genome_interval_tree.addi(contig, start, end, 
                                               data=AnnotatedInterval(contig, start, end, annotations))
            
    def get_overlaps(self, interval):
        tree_query = self.genome_interval_tree[interval.contig].search(interval.start, interval.end)
        return [found_interval.data for found_interval in tree_query]
    
    def __iter__(self):
        return chain(*[self.genome_interval_tree[contig].__iter__() for contig in self.contig_dict.keys()])

    def __str__(self):
        return str(self.annotated_intervals_df)

    def __repr__(self):
        return self.__str__()
    
    def __len__(self):
        return len(self.annotated_intervals_df)


def plot_hist(x_y_weights, bins, xlabel, ylabel, cutoff, vmin=10**0, vmax=10**8, title=None, output_file=None):
    x_y_weights_ary = np.array(x_y_weights)
    x = x_y_weights_ary[:, 0]
    y = x_y_weights_ary[:, 1]
    weights = x_y_weights_ary[:, 2]
    
    nullfmt = NullFormatter()         # no tick labels
    
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_hist2d = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.3]
    rect_histy = [left_h, bottom, 0.3, height]

    # start with a rectangular Figure
    plt.figure(figsize=(8, 8))

    axHist2d = plt.axes(rect_hist2d)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)
    
    axHistx.set_title(title, fontsize=16)

    # no tick labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # the scatter plot:
    fig = axHist2d.hist2d(x, y, 
                          weights=weights,
                          bins=[bins, bins],
                          norm=LogNorm(),
                          vmin=vmin,
                          vmax=vmax,
                          cmap='YlGnBu')
    axHist2d.set_xlabel(xlabel, fontsize=16)
    axHist2d.set_ylabel(ylabel, fontsize=16)
    cbar = plt.colorbar(fig[3])
    cbar.set_label('Number of bases', fontsize=12, rotation=270, labelpad=10)

    axHist2d.set_xlim((bins[0], bins[-1]))
    axHist2d.set_ylim((bins[0], bins[-1]))
    axHist2d.plot(axHist2d.get_xlim(), axHist2d.get_ylim(), 
                  ls='dashed')

    axHist2d.plot([bins[0], bins[-1]-cutoff],
                  [bins[0]+cutoff, bins[-1]],
                  ls=":", color="lightgray")
    axHist2d.plot([bins[0]+cutoff, bins[-1]],
                  [bins[0], bins[-1]-cutoff],
                  ls=":", color="lightgray")

    axHistx.hist(x, 
                 weights=weights, 
                 bins=bins,
                 histtype='bar',
                 alpha=0.5,
                 log=True)
    axHisty.hist(y,
                 weights=weights, 
                 bins=bins,
                 histtype='bar',
                 alpha=0.5,
                 log=True,
                 orientation='horizontal')

    axHistx.set_xlim(axHist2d.get_xlim())
    axHistx.set_ylim([vmin, vmax])
    axHistx.set_ylabel('Number of bases', fontsize=12)
    axHisty.set_ylim(axHist2d.get_ylim())
    axHisty.set_xlim([vmin, vmax])
    axHisty.set_xlabel('Number of bases', fontsize=12)
    
    plt.savefig(output_file, bbox_inches='tight')


def main(prefix, input1_segs, input2_segs, input1_format, input2_format,
         cutoff, xlab, ylab, title, print_concordance, vmax_factor, discordant_coords):
    # Check SEG format
    seg_format = {"model_seg":["CONTIG", "START", "END", "LOG2_COPY_RATIO_POSTERIOR_50"],
                  "regular_seg":["Chromosome", "Start", "End", "Segment_Mean"]}
    if input1_format not in seg_format:
        sys.exit("--seg1-type: specify either 'regular_seg' or 'model_seg'")
    if input2_format not in seg_format:
        sys.exit("--seg2-type: specify either 'regular_seg' or 'model_seg'")

    fo = None
    if discordant_coords:
        fo = open(discordant_coords, "w")
        fo.write("chrom\tstart\tend\tabs_diff\tseg1_cr\tseg2_cr\n")

    cr_bins = np.arange(0, 2.0, 0.05) - 0.025
    length_threshold = 0
    cr_threshold = cutoff
    
    # Load seg files
    cnv_profile = []
    cr_concordant = 0
    total_overlap = 0
    for fi in range(len(input1_segs)):
        input1_aic = AnnotatedIntervalCollection(input1_segs[fi],
                                                 contig_header=seg_format[input1_format][0],
                                                 start_header=seg_format[input1_format][1],
                                                 end_header=seg_format[input1_format][2])
        input2_aic = AnnotatedIntervalCollection(input2_segs[fi], 
                                                 contig_header=seg_format[input2_format][0],
                                                 start_header=seg_format[input2_format][1],
                                                 end_header=seg_format[input2_format][2])
        for input1_segment in input1_aic:
            input1_cr = 2**input1_segment.data.annotations[seg_format[input1_format][3]]

            segments = input2_aic.get_overlaps(input1_segment.data)

            if len(segments) > 0:
                for segment in segments:
                    input2_cr = 2**segment.annotations[seg_format[input2_format][3]]
                    overlap = input1_segment.data.overlap_in_bp(segment)

                    cnv_profile.append([input1_cr, input2_cr, overlap])

                    if abs(input1_cr - input2_cr) < cr_threshold:
                        cr_concordant += overlap
                    else:
                        if discordant_coords:
                            fo.write("\t".join(map(str, input1_segment.data.overlap_coord(segment) +
                                                   [abs(input1_cr-input2_cr), input1_cr, input2_cr])) + "\n")
                    total_overlap += overlap
    plot_hist(cnv_profile,
              xlabel=xlab, ylabel=ylab, cutoff=cutoff,
              bins=cr_bins, vmin=1000, vmax=10**vmax_factor,
              title='{0}\nconcordant = {2}\nbases = {1} / {3}'.format(title,
                    cr_concordant,
                    cr_concordant/float(total_overlap),
                    total_overlap),
              output_file='{0}-diff-cutoff-{1}.svg'.format(prefix, cutoff))
    if discordant_coords:
        fo.close()

    if print_concordance:
        input1_sample_name = os.path.basename(input1_segs[0])
        input2_sample_name = os.path.basename(input2_segs[0])
        print("seg1\tseg2\tbases_concordant_cnv\ttotal_bases_examined")
        print("\t".join([input1_sample_name, input2_sample_name,
                         str(cr_concordant), str(total_overlap),
                         str(cr_concordant/float(total_overlap))]))
    


if __name__ == "__main__":
    parser = ArgumentParser(description="Check CNV concordance between SEGs")
    parser.add_argument("--vmax-factor", type=int, help="factor for number of bases - i.e. 10**x (default: 9)", default=9)
    parser.add_argument("--discordant-coords", help="Generate a file that contains discordant overlapping CNVs")
    parser.add_argument("--print-concordance", help="Print concordance level", action="store_true")
    parser.add_argument("--cutoff", type=float, help="cutoff for neutral and CNA in copy number space (default: 0.1)", default=0.1)
    parser.add_argument("--seg1-type", help="group 1 tumor seg file format [regular_seg / model_seg]", required=True)
    parser.add_argument("--seg2-type", help="group 2 tumor seg file format [regular_seg / model_seg]", required=True)
    parser.add_argument("--seg1", nargs="+", help="group 1 tumor seg file(s)", required=True)
    parser.add_argument("--seg2", nargs="+", help="group 2 tumor seg file(s)", required=True)
    parser.add_argument("--prefix", help="output prefix", required=True)
    parser.add_argument("--title", help="plot title", default="-")
    parser.add_argument("--xlab", help="X axis label (seg1 name)", default="Input 1 copy ratio")
    parser.add_argument("--ylab", help="Y axis label (seg2 name)", default="Input 2 copy ratio")

    args = parser.parse_args()
    
    main(args.prefix, args.seg1, args.seg2, args.seg1_type, args.seg2_type,
         args.cutoff, args.xlab, args.ylab, args.title, args.print_concordance,
         args.vmax_factor, args.discordant_coords)
