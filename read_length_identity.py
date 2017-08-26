#!/usr/bin/env python3
"""
This script produces a table with information for each read using a minimap2 PAF as input:
  * length
  * identity
  * relative length

If less than half of a read aligns, it is deemed unaligned and given an identity of 0. If more than
half aligns, only the aligned parts are used to determined the read identity.

Relative length is included to see if the reads are systematically too short or too long.
"""


import sys
import collections
import statistics


def main():
    read_lengths = {}
    read_alignments = collections.defaultdict(list)

    paf_filename = sys.argv[1]
    with open(paf_filename, 'rt') as paf:
        for line in paf:
            paf_parts = line.strip().split('\t')
            if len(paf_parts) < 11:
                continue

            read_name = paf_parts[0]
            read_length = int(paf_parts[1])
            read_lengths[read_name] = read_length

            read_start = int(paf_parts[2])
            read_end = int(paf_parts[3])
            ref_start = int(paf_parts[7])
            ref_end = int(paf_parts[8])
            identity = 100.0 * int(paf_parts[9]) / int(paf_parts[10])

            read_alignments[read_name].append((read_start, read_end, ref_start, ref_end, identity))

    print('\t'.join(['Name', 'Length', 'Identity', 'Relative length']))
    for read_name, read_length in read_lengths.items():
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

        relative_length = 100.0 * total_read_length / total_ref_length

        print('\t'.join([read_name, str(read_length), str(whole_read_identity), str(relative_length)]))


if __name__ == '__main__':
    main()
