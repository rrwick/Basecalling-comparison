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

This script takes an assembly as input and produces an output of 'reads': the assembly chopped into
overlapping pieces. It is used in assembly polishing to turn individual Racon-polished assemblies
into 'reads' which can then be used as input for a final Racon polish.
"""

import random
import sys
import uuid


def main():
    assembly = load_fasta(sys.argv[1])
    random.seed(int(sys.argv[2]))
    read_length = int(sys.argv[3])
    assert len(assembly) == 1
    assembly_seq = assembly[0][1]
    assert len(assembly_seq) > 1000000

    # Add a bit of overlap so reads can span the junction.
    assembly_seq = assembly_seq + assembly_seq[:read_length]

    read_start = 0
    read_end = read_length
    while read_end <= len(assembly_seq):
        read_seq = assembly_seq[read_start:read_end]
        if random.random() < 0.5:
            read_seq = reverse_complement(read_seq)
        read_qual = ''.join([chr(random.randint(65, 75)) for _ in range(read_length)])
        print('@' + str(uuid.uuid4()))
        print(read_seq)
        print('+')
        print(read_qual)
        read_start += read_length // 3
        read_end += read_length // 3


def load_fasta(fasta_filename):
    fasta_seqs = []
    with open(fasta_filename, 'rt') as fasta_file:
        name = ''
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    fasta_seqs.append((name.split()[0], ''.join(sequence), name))
                    sequence = []
                name = line[1:]
            else:
                sequence.append(line)
        if name:
            fasta_seqs.append((name.split()[0], ''.join(sequence), name))
    return fasta_seqs


REV_COMP_DICT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                 'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
                 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
                 'D': 'H', 'H': 'D', 'N': 'N',
                 'r': 'y', 'y': 'r', 's': 's', 'w': 'w',
                 'k': 'm', 'm': 'k', 'b': 'v', 'v': 'b',
                 'd': 'h', 'h': 'd', 'n': 'n',
                 '.': '.', '-': '-', '?': '?'}


def complement_base(base):
    try:
        return REV_COMP_DICT[base]
    except KeyError:
        return 'N'


def reverse_complement(seq):
    return ''.join([complement_base(x) for x in seq][::-1])


if __name__ == '__main__':
    main()
