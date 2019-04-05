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

# This script will execute each version of Flappie on the reads.


# Edit the following paths before running, as appropriate for your environment:
fast5_dir=/path/to/fast5s                  # all fast5s must be directly in this directory (i.e. not nested in subdirectories)
output_prefix=/path/to/basecalling_output  # Flappie output fastqs will go in this directory
binaries=/path/to/binaries                 # should contain Flappie binaries named with the version number, e.g. flappie-v1.1


cd $output_prefix
export OPENBLAS_NUM_THREADS=1

for v in 1.0 1.1; do
    out_dir="$output_prefix"/flappie_v"$v"
    flappie="$binaries"/flappie-"$v"

    # I can't use time for Flappie because I run it in separate processes with the parallel
    # command, so I just run date before and after to get the total time.
    date
    find $fast5_dir -name \*.fast5 | parallel -P 12 -X $flappie --model r941_native > flappie_v"$v".fastq
    date
done
