#!/usr/bin/env python3
"""
Copyright 2019 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Basecalling-comparison

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should
have received a copy of the GNU General Public License along with this program. If not, see
<http://www.gnu.org/licenses/>.

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
