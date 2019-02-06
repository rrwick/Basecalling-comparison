#!/usr/bin/env python3
"""
This script takes an assembly as input and produces an output of 'reads': the assembly chopped into pieces.
It is for assessing the distribution of identity over the assembly.
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
