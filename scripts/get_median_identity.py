#!/usr/bin/env python3
"""
This script reads a tsv file (for either read or assembly data) and returns the median identity.
It can optionally take a total number of sequences, which will make it add zeros to get to that
number (useful for basecallers that didn't basecall all of the reads).
"""

import statistics
import sys

identities = []
with open(sys.argv[1], 'rt') as data_file:
    for line in data_file:
        parts = line.strip().split()
        if parts[0] == 'Name':
            continue
        identities.append(float(parts[2]))

try:
    total_count = int(sys.argv[2])
    while len(identities) < total_count:
        identities.append(0.0)
except IndexError:
    pass

print(statistics.median(identities))
