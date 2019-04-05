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

This script filters fast5 reads based on their alignment to reference contigs.
The reads which pass filtering should be suitable for Sloika training.
"""

import argparse
import glob
import os
import re
import shutil
import sys


def get_arguments():
    parser = argparse.ArgumentParser(description='Trim fast5 files at the signal level')

    parser.add_argument('--min_basecalled_length', type=int, required=False, default=5000,
                        help='Reads will fail if basecalled to shorter than this amount')
    parser.add_argument('--max_unaligned_bases', type=int, required=False, default=100,
                        help='Reads will fail if they have more than this many bases unaligned')
    parser.add_argument('--max_window_indels', type=float, required=False, default=0.8,
                        help='Reads will fail if they have any CIGAR window with more than this '
                             'fraction of indels')
    parser.add_argument('--window_size', type=int, required=False, default=25,
                        help='Window size for CIGAR indel check')

    parser.add_argument('in_fast5_dir', type=str,
                        help='Directory containing trimmed fast5 files')
    parser.add_argument('seq_summary', type=str,
                        help='Albacore/Guppy\'s sequencing_summary.txt file for the basecalling '
                             'of the trimmed fast5 files')
    parser.add_argument('reference', type=str,
                        help='Reference FASTA (used to generate alignment)')
    parser.add_argument('paf_alignment', type=str,
                        help='Minimap2 alignment of basecalled reads to reference contigs')
    parser.add_argument('out_fast5_dir', type=str,
                        help='Output directory of passing fast5 files')
    parser.add_argument('out_ref', type=str,
                        help='Output file per-read reference sequences')

    args = parser.parse_args()
    return args


def main():
    args = get_arguments()
    filename_to_read_id = read_seq_summary(args.seq_summary)
    read_id_to_alignment_info = read_paf(args.paf_alignment)
    make_output_dir(args.out_fast5_dir)
    references = load_fasta(args.reference)

    with open(args.out_ref, 'w') as out_ref:
        for fast5_file in glob.iglob(args.in_fast5_dir + '/**/*.fast5', recursive=True):
            print()
            print(fast5_file)
            filename = fast5_file.rpartition('/')[2]
            read_id = filename_to_read_id[filename]
            print('    ID: {}'.format(read_id))

            try:
                read_length, read_start, read_end, cigar, _, strand, contig_name, contig_start, contig_end = read_id_to_alignment_info[read_id]
            except KeyError:
                print('    FAIL due to no alignment')
                continue

            print('    Length: {:,} bp'.format(read_length))
            if read_length < args.min_basecalled_length:
                print('    FAIL due to short length')
                continue

            unaligned_bases = read_start + (read_length - read_end)
            print('    Unaligned: {:,} bp'.format(unaligned_bases))
            if unaligned_bases > args.max_unaligned_bases:
                print('    FAIL due to too much unaligned')
                continue

            expanded_cigar = get_expanded_cigar(cigar)
            worst_window_indel_fraction = 0.0
            for i in range(len(expanded_cigar) - args.window_size):
                cigar_window = expanded_cigar[i:i+args.window_size]
                window_indel_fraction = (cigar_window.count('I') + cigar_window.count('D')) / args.window_size
                worst_window_indel_fraction = max(worst_window_indel_fraction, window_indel_fraction)
            print('    Worst window indels: {:.4f}'.format(worst_window_indel_fraction))
            if worst_window_indel_fraction > args.max_window_indels:
                print('    FAIL due to bad window indels')
                continue

            print('    PASS')
            new_file = '{}/{}'.format(args.out_fast5_dir, filename)
            print('    making a copy here: {}'.format(new_file))
            if os.path.exists(new_file):
                print()
                sys.exit('Error: {} already exists'.format(new_file))
            shutil.copyfile(fast5_file, new_file)

            print('    getting reference sequence from {} ({}-{}, {} strand)'.format(contig_name, contig_start, contig_end, strand))
            ref_header = '>' + os.path.basename(os.path.splitext(fast5_file)[0])
            ref_seq = references[contig_name][contig_start:contig_end]
            if strand == '-':
                ref_seq = reverse_complement(ref_seq)

            out_ref.write(ref_header)
            out_ref.write('\n')
            out_ref.write(ref_seq)
            out_ref.write('\n')


def read_seq_summary(seq_summary_filename):
    filename_to_read_id = {}
    with open(seq_summary_filename, 'r') as seq_summary:
        for line in seq_summary:
            parts = line.split('\t')
            if parts[0] == 'filename':  # header line
                continue
            filename, read_id = parts[0], parts[1]
            filename_to_read_id[filename] = read_id
    return filename_to_read_id


def read_paf(paf_filename):
    read_id_to_alignment_info = {}
    with open(paf_filename, 'r') as paf_file:
        for line in paf_file:
            parts = line.split('\t')
            read_id = parts[0]
            read_length = int(parts[1])
            read_start = int(parts[2])
            read_end = int(parts[3])
            strand = parts[4]
            contig_name = parts[5]
            contig_start = int(parts[7])
            contig_end = int(parts[8])
            cigar, alignment_score = None, None
            for part in parts:
                if part.startswith('cg:Z:'):
                    cigar = part[5:]
                if part.startswith('AS:i:'):
                    alignment_score = int(part[5:])
            assert cigar is not None and alignment_score is not None
            alignment_info = read_length, read_start, read_end, cigar, alignment_score, strand, contig_name, contig_start, contig_end

            # If there are multiple alignments for a read, we keep the one with the best score.
            if read_id in read_id_to_alignment_info:
                _, _, _, _, existing_score, _, _, _, _ = read_id_to_alignment_info[read_id]
                if alignment_score > existing_score:
                    read_id_to_alignment_info[read_id] = alignment_info
            else:
                read_id_to_alignment_info[read_id] = alignment_info
    return read_id_to_alignment_info


def get_expanded_cigar(cigar):
    expanded_cigar = []
    cigar_parts = re.findall(r'\d+\w', cigar)
    for cigar_part in cigar_parts:
        num = int(cigar_part[:-1])
        letter = cigar_part[-1]
        if letter == 'M' or letter == 'I' or letter == 'D':
            expanded_cigar.append(letter * num)
    return ''.join(expanded_cigar)


def make_output_dir(out_dir):
    if os.path.exists(out_dir):
        sys.exit('Error: output directory already exists')
    os.makedirs(out_dir)


def load_fasta(fasta_filename):
    fasta_seqs = {}
    with open(fasta_filename, 'rt') as fasta_file:
        name = ''
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    fasta_seqs[name.split()[0]] = ''.join(sequence)
                    sequence = []
                name = line[1:]
            else:
                sequence.append(line)
        if name:
            fasta_seqs[name.split()[0]] = ''.join(sequence)
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
