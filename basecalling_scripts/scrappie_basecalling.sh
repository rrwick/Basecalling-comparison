#!/usr/bin/env bash

# Copyright 2019 Ryan Wick (rrwick@gmail.com)
# https://github.com/rrwick/Basecalling-comparison

# This program is free software: you can redistribute it and/or modify it under the terms of the
# GNU General Public License as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version. This program is distributed in the hope that it
# will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You
# should have received a copy of the GNU General Public License along with this program. If not,
# see <http://www.gnu.org/licenses/>.

# This script will execute each version of Scrappie on the reads.


# Edit the following paths before running, as appropriate for your environment:
fast5_dir=/path/to/fast5s                  # all fast5s must be directly in this directory (i.e. not nested in subdirectories)
output_prefix=/path/to/basecalling_output  # Scrappie output fastas will go in this directory
binaries=/path/to/binaries                 # should contain Scrappie binaries named with the version number, e.g. scrappie-v1.4.1


cd $output_prefix
export OMP_NUM_THREADS=$(nproc)
export OPENBLAS_NUM_THREADS=1

# For recent versions of Scrappie, I ran both raw and event basecalling:
for v in 1.4.1 1.4.0 1.3.2 1.3.1 1.3.0 1.2.0 1.1.1 1.1.0; do
    out_dir="$output_prefix"/scrappie_v"$v"
    scrappie="$binaries"/scrappie-v"$v"
    /usr/bin/time -v -o scrappie_raw_v"$v".time $scrappie raw $fast5_dir > scrappie_raw_v"$v".fasta
    /usr/bin/time -v -o scrappie_events_v"$v".time $scrappie events $fast5_dir > scrappie_events_v"$v".fasta
done

# These versions can't call events on their own - they rely on events already
# being present, e.g. from Albacore. I therefore skipped scrappie events
# (because it's not a stand-alone basecaller) and just ran scrappie raw.
for v in 1.0.0 0.3.2 0.3.1 0.3.0; do
    out_dir="$output_prefix"/scrappie_v"$v"
    scrappie="$binaries"/scrappie-v"$v"
    /usr/bin/time -v -o scrappie_raw_v"$v".time $scrappie raw $fast5_dir > scrappie_raw_v"$v".fasta
done


# Versions earlier than 0.3.0 can only do event-based basecalling, so I'm skipping all of them.
