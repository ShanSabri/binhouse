#!/usr/bin/env python

'''
bed_to_fpkm.py
Shan Sabri
Dec. 14, 2018

Compute FPKM values from a standard BED file
I: BED file
O: Wiggle file
'''

import gzip
from optparse import OptionParser
import os
import sys
import datetime
import math


def mean(array):
    return sum(array) / float(len(array))

def std(array):
    m = mean(array)
    return math.sqrt(sum((x - m) ** 2 for x in array) / float(len(array)))

def echo(message):
    print "[%s] %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), str(message))

def overlap(s1, e1, s2, e2):
    return max(min(e1, e2) - max(s1, s2), 0)

def get_fpkms(bed_fname,
              chrom_sizes,
              filter_chrom=None,
              bin_size=200,
              extend=200,
              shift=0,
              min_fpkm=0,
              normalize_factor=None,
              var_normalize=False):

    if filter_chrom and type(filter_chrom) is not list:
        filter_chrom = [filter_chrom]

    # fpkms = dict((chrom, [0] * (1 + (chrom_sizes[chrom] - 1) / BIN_SIZE)) for chrom in (filter_chrom or chrom_sizes))
    fpkms = dict((chrom, [0] * (chrom_sizes[chrom] / bin_size)) for chrom in (filter_chrom or chrom_sizes))

    total_reads = 0
    chrom_total = dict((chrom, 0) for chrom in fpkms)

    echo('Counting reads in ' + bed_fname)
    with gzip.open(bed_fname) if bed_fname.endswith('.gz') else open(bed_fname) as bed_f:
        for line in bed_f:
            buf = line.strip().split()
            if len(buf) == 4:
                chrom, start, end, strand = buf
            else:
                chrom, start, end = buf[:3]
                strand = buf[5]

            if filter_chrom and chrom not in filter_chrom:
                continue

            start = int(start)
            end = int(end)

            if strand == '+':
                start += shift

                if extend > 0:
                    end = min(start + extend, chrom_sizes[chrom])
            else:

                end -= shift

                if extend > 0:
                    start = max(0, end - extend)

            effective_read_length = end - start

            start_bin = start / bin_size
            end_bin = min(1 + end / bin_size, len(fpkms[chrom]))

            for bin_idx in xrange(start_bin, end_bin):
                fraction = float(overlap(start, end,
                                         bin_idx * bin_size, (bin_idx + 1) * bin_size)) / effective_read_length
                # fraction = 1
                fpkms[chrom][bin_idx] += fraction

            total_reads += 1
            chrom_total[chrom] += 1

    factor = float(10 ** 9 / bin_size) / total_reads

    if normalize_factor is not None:
        echo('External normalizing factor: ' + str(normalize_factor))
        factor *= normalize_factor

    if var_normalize:
        stddev = std([n_reads * factor if n_reads * factor >= min_fpkm else 0
                          for chrom in fpkms if chrom_total[chrom] > 0
                            for n_reads in fpkms[chrom]])
        echo('Stddev: ' + str(stddev) + '. Normalizing to variance 1!')
        factor /= stddev

    return dict((chrom, [n_reads * factor if n_reads * factor >= min_fpkm else 0
                         for n_reads in fpkms[chrom]])
                for chrom in fpkms if chrom_total[chrom] > 0)


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-i", "--input", type="string", dest="in_fname",
                      help="Input bed file with aligned reads (can be gzipped bed)", metavar="FILE")
    parser.add_option("-g", "--genome", type="string", dest="genome_fname",
                      help="Genome file with chromosome lengths [ex. mm9.txt]", metavar="FILE")

    parser.add_option("-o", "--output", type="string", dest="out_fname",
                      help="Output file name. The output is in wiggle format.", metavar="FILE")

    # parser.add_option("-z", "--gzip", action="store_true", dest="gzip", default=False,
    #                   help="Gzip the output. [%default]")

    parser.add_option("--variable-step", action="store_true", dest="var_step", default=False,
                      help="Output variable step wiggle format. [%default]")

    parser.add_option("-b", "--bin-size", type="int", dest="bin_size", default=200,
                      help="Bin size in base pairs [%default]", metavar="INT")

    parser.add_option("-e", "--extend", type="int", dest="extend", default=0,
                      help="Extend the reads so many base pairs from the start of the alignment [%default]",
                      metavar="INT")

    parser.add_option("-s", "--shift", type="int", dest="shift", default=0,
                      help="Shifts the start of each read by this many base pairs in the direction of alignment. "
                           "[%default]",
                      metavar="INT")

    parser.add_option("-m", "--min-fpkm", type="int", dest="min_fpkm", default=0,
                      help="FPKM values less than that will be set to 0. [%default]",
                      metavar="INT")

    parser.add_option("--variance-normalize", action="store_true", dest="var_normalize", default=False,
                      help="Normalize to variance 1. [%default]")

    parser.add_option("--normalize-factor", type="float", dest="normalize_factor", default=None,
                      help="External normalizing factor that will be multiplied in with each value. [%default]",
                      metavar="FLOAT")

    parser.add_option("--smooth", type="int", dest="smooth", default=0,
                      help="Smooth the values with +- this amount of bins. [%default]",
                      metavar="INT")

    parser.add_option("-c", "--chromosomes", type="string", dest="filter_chrom", default='',
                      help="List of chromosomes to output separated by commas (ex. chr1,chr2,chr3). "
                           "Default: No filtering",
                      metavar="STRING")

    (options, args) = parser.parse_args()

    # if no options were given by the user, print help and exit
    if len(sys.argv) == 1:
        parser.print_help()
        exit(0)

    # calculate genome-wide FPKM values

    filter_chrom = None if not options.filter_chrom else options.filter_chrom.split(',')
    echo('Filtering: ' + str(filter_chrom))

    with open(options.genome_fname) as genome_f:
        chrom_sizes = dict((l[0], int(l[1])) for l in map(lambda x: x.strip().split(), genome_f))

    fpkms = get_fpkms(options.in_fname,
                      chrom_sizes,
                      filter_chrom,
                      options.bin_size,
                      options.extend,
                      options.shift,
                      options.min_fpkm,
                      options.normalize_factor,
                      options.var_normalize
                      )

    echo('Storing output in ' + options.out_fname)

    with (gzip.open(options.out_fname, 'w') if options.out_fname.endswith('.gz') else open(options.out_fname, 'w')) as out_f:
        title = os.path.split(options.out_fname)[1].replace('.bed', '').replace('.gz', '').replace('.wig', '')

        out_f.write('track type=wiggle_0 name="%s" description="%s"\n' % (title, title))

        for chrom in sorted(fpkms):
            out_f.write('%s\tchrom=%s\tstart=1\tstep=%d\tspan=%d\n' % ('variableStep' if options.var_step else 'fixedStep',
                                                                       chrom,
                                                                       options.bin_size,
                                                                       options.bin_size))
            for bin_idx in xrange(len(fpkms[chrom])):
                out_f.write(('%d\t' % (bin_idx * options.bin_size) if options.var_step else '') +
                            '%.2lf\n' % mean(fpkms[chrom][max(0, bin_idx - options.smooth):
                                                          min(len(fpkms[chrom]), bin_idx + options.smooth + 1)]))

    echo('Done')
