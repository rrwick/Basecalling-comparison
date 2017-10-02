#!/usr/bin/env python3
"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Basecalling-comparison

This script produces a table with information for each read.

Inputs:
  * Read file
  * minimap2 PAF

Output:
  * tsv file with these columns: length, identity, relative length

If less than half of a read aligns, it is deemed unaligned and given an identity of 0. If more than
half aligns, only the aligned parts are used to determined the read identity.

Relative length is included to see if the reads are systematically too short or too long.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should
have received a copy of the GNU General Public License along with this program. If not, see
<http://www.gnu.org/licenses/>.
"""


import sys
import collections
import statistics
import os
import gzip


def main():
    read_alignments = collections.defaultdict(list)

    read_filename = sys.argv[1]
    paf_filename = sys.argv[2]

    read_lengths = get_read_lengths(read_filename)

    with open(paf_filename, 'rt') as paf:
        for line in paf:
            paf_parts = line.strip().split('\t')
            if len(paf_parts) < 11:
                continue

            read_name = paf_parts[0]
            read_length = int(paf_parts[1])
            assert read_length == read_lengths[read_name]

            read_start = int(paf_parts[2])
            read_end = int(paf_parts[3])
            ref_start = int(paf_parts[7])
            ref_end = int(paf_parts[8])
            identity = 100.0 * int(paf_parts[9]) / int(paf_parts[10])

            read_alignments[read_name].append((read_start, read_end, ref_start, ref_end, identity))

    print('\t'.join(['Name', 'Length', 'Identity', 'Relative length']))
    read_names = sorted(read_lengths.keys())
    for read_name in read_names:
        read_length = read_lengths[read_name]
        alignments = read_alignments[read_name]
        identity_by_base = [0.0] * read_length

        total_read_length = 0
        total_ref_length = 0

        for read_start, read_end, ref_start, ref_end, identity in alignments:
            for i in range(read_start, read_end):
                if identity > identity_by_base[i]:
                    identity_by_base[i] = identity
            total_read_length += read_end - read_start
            total_ref_length += ref_end - ref_start

        # If less than half the read aligned, then we call it an unaligned read.
        if identity_by_base.count(0.0) > read_length / 2.0:
            whole_read_identity = 0.0

        # Otherwise, the read's identity is the average of the aligned parts.
        else:
            whole_read_identity = statistics.mean([x for x in identity_by_base if x > 0.0])

        if whole_read_identity > 0.0:
            relative_length = str(100.0 * total_read_length / total_ref_length)
        else:
            relative_length = ''

        print('\t'.join([read_name, str(read_length), str(whole_read_identity), relative_length]))


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for filetype, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = filetype
    if compression_type == 'bz2':
        sys.exit('Error: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('Error: cannot use zip format - use gzip instead')
    return compression_type


def get_sequence_file_type(filename):
    """
    Determines whether a file is FASTA or FASTQ.
    """
    if not os.path.isfile(filename):
        sys.exit('Error: could not find ' + filename)
    if get_compression_type(filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open

    with open_func(filename, 'rt') as seq_file:
        try:
            first_char = seq_file.read(1)
        except UnicodeDecodeError:
            first_char = ''

    if first_char == '>':
        return 'FASTA'
    elif first_char == '@':
        return 'FASTQ'
    else:
        raise ValueError('File is neither FASTA or FASTQ')


def get_read_lengths(filename):
    """
    Returns a dictionary of read names to read lengths.
    """
    try:
        file_type = get_sequence_file_type(filename)
        if file_type == 'FASTA':
            return get_fasta_lengths(filename)
        else:  # FASTQ
            return get_fastq_lengths(filename)
    except IndexError:
        sys.exit('\nError: ' + filename + ' could not be parsed - is it formatted correctly?')


def get_fasta_lengths(fasta_filename):
    if get_compression_type(fasta_filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    read_lengths = {}
    with open_func(fasta_filename, 'rt') as fasta_file:
        name = ''
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    read_lengths[name.split()[0]] = len(''.join(sequence))
                    sequence = []
                name = line[1:]
            else:
                sequence.append(line)
        if name:
            read_lengths[name.split()[0]] = len(''.join(sequence))
    return read_lengths


def get_fastq_lengths(fastq_filename):
    """
    Returns a list of tuples (header, seq) for each record in the fastq file.
    """
    if get_compression_type(fastq_filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    read_lengths = {}
    with open_func(fastq_filename, 'rt') as fastq:
        for line in fastq:
            name = line.strip()[1:].split()[0]
            sequence = next(fastq).strip()
            next(fastq)
            next(fastq)
            read_lengths[name] = len(sequence)
    return read_lengths


if __name__ == '__main__':
    main()
