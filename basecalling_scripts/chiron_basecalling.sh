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

# This script will execute each version of Chiron on the reads. Chiron takes a long time to run,
# so I ran it on smaller batches of files so I could start/stop the basecalling as needed.


# Edit the following paths before running, as appropriate for your environment:
fast5_dir=/path/to/fast5s_100_read_dirs    # should contain numbered subdirectories (000 to 151), each containing fast5 files
output_prefix=/path/to/basecalling_output  # Chiron output fastqs will go in this directory
chiron_dirs=/path/to/Chiron                # should contain a subdirectory for each version of Chiron, e.g. Chiron-v0.4.2

# TensorFlow/CUDA requirements might change between Chiron versions:
cuda_bin_path=/path/to/cuda-8.0/bin
cuda_lib_path=/path/to/cuda-8.0/lib64


cd $output_prefix

# Versions 0.4, 0.4.1 and 0.4.2 all seem to produce the same reads, so I'm only running v0.4.2.
v=0.4.2
for d in {000..151}; do
    out_dir="$output_prefix"/chiron_v"$v"_"$d"
    out_file=chiron_v"$v"_"$d".fastq
    time_file=chiron_v"$v"_"$d".time
    if [ ! -f $out_file ]; then
        rm -rf $out_dir
        rm -f $time_file
        /usr/bin/time -v -o $time_file python "$chiron_dirs"/Chiron-"$v"/chiron/entry.py call -i "$fast5_dir"/"$d" -o $out_dir || break
        cat "$out_dir"/result/*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > $out_file
    fi
done
cat chiron_v"$v"_*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > chiron_v"$v".fastq
cat chiron_v"$v"_*.time > chiron_v"$v".time

# Older versions need TensorFlow 1.0.1, which needs CUDA 8.0.
export PATH="$cuda_bin_path":"$PATH"
export LD_LIBRARY_PATH="$cuda_lib_path":"$LD_LIBRARY_PATH"

for v in 0.3 0.2; do
    for d in "{000..151}"; do
        out_dir="$output_prefix"/chiron_v"$v"_"$d"
        out_file=chiron_v"$v"_"$d".fastq
        time_file=chiron_v"$v"_"$d".time
        if [ ! -f $out_file ]; then
            rm -rf $out_dir
            rm -f $time_file
            /usr/bin/time -v -o $time_file python "$chiron_dirs"/Chiron-"$v"/chiron/entry.py call -i "$fast5_dir"/"$d" -o $out_dir || break
            cat "$out_dir"/result/*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > $out_file
        fi
    done
    cat chiron_v"$v"_*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > chiron_v"$v".fastq
    cat chiron_v"$v"_*.time > chiron_v"$v".time
done
