#!/usr/bin/env python3
"""
This script takes as input the table from the read_length_identity.py script. It produces two
tables for making histograms: one for read identity and one for relative length.
"""

import sys
import collections


def main():
    input_table_filename = sys.argv[1]
    identity_bin_size = float(sys.argv[2])
    relative_length_bin_size = float(sys.argv[3])
    identity_histogram_filename = sys.argv[4]
    relative_length_histogram_filename = sys.argv[5]

    identity_bins = collections.defaultdict(int)
    relative_length_bins = collections.defaultdict(int)

    with open(input_table_filename, 'rt') as input_table:
        for line in input_table:
            line_parts = line.rstrip().split()
            if line_parts[0] == 'Name':
                continue
            name = line_parts[0]
            length = int(line_parts[1])
            identity = float(line_parts[2])
            relative_length = float(line_parts[3])

            identity_bins[round_down(identity, identity_bin_size)] += length
            relative_length_bins[round_down(relative_length, relative_length_bin_size)] += length

    write_histogram_to_file(identity_histogram_filename, identity_bins, identity_bin_size)
    write_histogram_to_file(relative_length_histogram_filename, relative_length_bins, relative_length_bin_size)


def write_histogram_to_file(filename, bins, bin_size):
    min_bin = min(bins.keys())
    max_bin = max(bins.keys())
    with open(filename, 'wt') as histogram_file:
        histogram_file.write('Bin\tValue\n')
        for i in range(min_bin, max_bin + bin_size):
            histogram_file.write(str(i))
            histogram_file.write('\t')
            histogram_file.write(str(bins[i]))
            histogram_file.write('\n')


def round_down(num, step):
    result = int(num / step) * step
    return float('%.4f' % result)  # deal with floating point error


if __name__ == '__main__':
    main()
