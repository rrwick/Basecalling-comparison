#!/usr/bin/env python3
"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Basecalling-comparison

This script adjusts read headers to be consistent between basecallers and compatible with
Nanopolish. After running, each read header should be in this format:
5a8d447e-84e2-4f6f-922c-5ad7269f688c_Basecall_1D_template 5210_N125509_20170425_FN2002039725_MN19691_sequencing_run_klebs_033_restart_87298_ch152_read14914_strand

It also sorts the reads alphabetically by their new headers and removes 0-length reads.

Usage:
fix_read_names.py input_reads.fastq.gz read_id_to_fast5 | gzip > output_reads.fastq.gz

It can take either fasta or fastq input and will output in the same format.

The read_id_to_fast5 file is a tab-delimited file with read IDs in the first column and fast5
filenames in the second. For example:
0000974e-e5b3-4fc2-8fa5-af721637e66c_Basecall_1D_template	5210_N125509_20170425_FN2002039725_MN19691_sequencing_run_klebs_033_restart_87298_ch173_read25236_strand.fast5
00019174-2937-4e85-b104-0e524d8a7ba7_Basecall_1D_template	5210_N125509_20170424_FN2002039725_MN19691_sequencing_run_klebs_033_75349_ch85_read2360_strand.fast5
000196f6-6041-49a5-9724-77e9d117edbe_Basecall_1D_template	5210_N125509_20170425_FN2002039725_MN19691_sequencing_run_klebs_033_restart_87298_ch200_read1975_strand.fast5

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should
have received a copy of the GNU General Public License along with this program. If not, see
<http://www.gnu.org/licenses/>.
"""


import sys
import os
import gzip
import re


def main():
    input_reads_filename = sys.argv[1]
    read_id_to_fast5_filename = sys.argv[2]

    print('\nReading ' + read_id_to_fast5_filename, file=sys.stderr, flush=True)
    read_id_to_fast5, fast5_to_read_id = {}, {}
    with open(read_id_to_fast5_filename, 'rt') as read_id_to_fast5_file:
        for line in read_id_to_fast5_file:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                read_id, fast5 = parts
                if fast5.endswith('.fast5'):
                    fast5 = fast5[:-6]
                if read_id in read_id_to_fast5:
                    sys.exit('Error: duplicate read ID in ' + read_id_to_fast5_filename + ': ' + read_id)
                if fast5 in fast5_to_read_id:
                    sys.exit('Error: duplicate fast5 in ' + read_id_to_fast5_filename + ': ' + fast5)
                read_id_to_fast5[read_id] = fast5
                fast5_to_read_id[fast5] = read_id

    print('Loading reads from ' + input_reads_filename, file=sys.stderr, flush=True)
    reads, read_type = load_fasta_or_fastq(input_reads_filename)

    print('Fixing read names', file=sys.stderr, flush=True)
    output_reads = []
    for header, seq, qual in reads:
        read_id, fast5_name = None, None
        try:
            read_id = re.search(r'[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}', header).group(0)
        except AttributeError:
            read_id = None
        try:
            fast5_name = re.search(r'\w+_ch\d+_read\d+_\w+', header).group(0)
        except AttributeError:
            fast5_name = None
        if read_id is None and fast5_name is None:
            sys.exit('Error: could not parse read header\n' + header)

        if read_id is not None:
            new_header = read_id + ' ' + read_id_to_fast5[read_id]
        else:
            new_header = fast5_to_read_id[fast5_name] + ' ' + fast5_name

        output_reads.append((new_header, seq, qual))

    print('Sorting reads', file=sys.stderr, flush=True)
    output_reads = sorted(output_reads)

    print('Outputting reads', file=sys.stderr, flush=True)
    for header, seq, qual in output_reads:
        if len(seq) == 0:
            continue
        if read_type == 'FASTA':
            print('>' + header)
            print(seq)
        else:  # read_type == 'FASTQ'
            print('@' + header)
            print(seq)
            print('+')
            print(qual)
    print('Done!\n', file=sys.stderr, flush=True)


def get_compression_type(filename):
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


def load_fasta_or_fastq(filename):
    try:
        file_type = get_sequence_file_type(filename)
        if file_type == 'FASTA':
            return load_fasta(filename), 'FASTA'
        else:  # FASTQ
            return load_fastq(filename), 'FASTQ'
    except IndexError:
        sys.exit('\nError: ' + filename + ' could not be parsed - is it formatted correctly?')


def load_fasta(fasta_filename):
    if get_compression_type(fasta_filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    fasta_seqs = []
    with open_func(fasta_filename, 'rt') as fasta_file:
        name = ''
        sequence = ''
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    fasta_seqs.append((name, sequence, ''))
                    sequence = ''
                name = line[1:]
            else:
                sequence += line
        if name:
            fasta_seqs.append((name, sequence, ''))
    return fasta_seqs


def load_fastq(fastq_filename):
    if get_compression_type(fastq_filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    reads = []

    extra_empty_line = False
    with open_func(fastq_filename, 'rt') as fastq:
        for line in fastq:
            _ = next(fastq)
            _ = next(fastq)
            _ = next(fastq)
            fifth_line = next(fastq).strip()
            if len(fifth_line) == 0:
                extra_empty_line = True
            break

    with open_func(fastq_filename, 'rt') as fastq:
        for line in fastq:
            name = line.strip()[1:]
            sequence = next(fastq).strip()
            _ = next(fastq)
            qualities = next(fastq).strip()
            if extra_empty_line:
                _ = next(fastq)
            reads.append((name, sequence, qualities))
    return reads


if __name__ == '__main__':
    main()
