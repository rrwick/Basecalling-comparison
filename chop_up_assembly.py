#!/usr/bin/env python3
"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Basecalling-comparison

This script takes an assembly as input and produces an output of 'reads': the assembly chopped into
pieces. It is for assessing the distribution of identity over the assembly.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should
have received a copy of the GNU General Public License along with this program. If not, see
<http://www.gnu.org/licenses/>.
"""

import sys


def main():
    assembly_filename = sys.argv[1]
    piece_size = int(sys.argv[2])

    contigs = load_fasta(assembly_filename)

    read_num = 0
    for contig in contigs:
        seq = contig[1]
        for i in range(0, len(seq), piece_size):
            piece_seq = seq[i:i+piece_size]
            if len(piece_seq) == piece_size:
                read_num += 1
                print('>' + str(read_num))
                print(piece_seq)


def load_fasta(fasta_filename):
    fasta_seqs = []
    with open(fasta_filename, 'rt') as fasta_file:
        name = ''
        sequence = ''
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    fasta_seqs.append((name.split()[0], sequence, name))
                    sequence = ''
                name = line[1:]
            else:
                sequence += line
        if name:
            fasta_seqs.append((name.split()[0], sequence, name))
    return fasta_seqs


if __name__ == '__main__':
    main()
