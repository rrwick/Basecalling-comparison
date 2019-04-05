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

# This script will execute each version of Guppy on the reads. In order to also use the
# custom-trained models, the model files(holtlab_kp_r9.4_r9.4.1_nov_2018.jsn and
# holtlab_kp_big_r9.4_r9.4.1_nov_2018.jsn) must be put into Guppy's data directory.


# Edit the following paths before running, as appropriate for your environment:
fast5_dir=/path/to/fast5s                   # all fast5s must be directly in this directory (i.e. not nested in subdirectories)
output_prefix=/path/to/basecalling_output   # Guppy output directories and consolidated fastqs will go in this directory
guppy_installers=/path/to/guppy_installers  # must contain the deb files for the Guppy versions - this script installs each version before running it


cd $output_prefix

# For the most recent versions of Guppy, I'm using both the default and flip-flop models:
for v in 2.1.3 2.2.3; do
    apt install "$guppy_installers"/ont-guppy_"$v"-1~xenial_amd64.deb
    out_dir="$output_prefix"/guppy_v"$v"
    /usr/bin/time -v -o guppy_v"$v".time guppy_basecaller --config dna_r9.4.1_450bps.cfg --device auto -i $fast5_dir -t 12 -s $out_dir
    cat "$out_dir"/*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > guppy_v"$v".fastq
    out_dir="$output_prefix"/guppy_v"$v"_flipflop
    /usr/bin/time -v -o guppy_v"$v"_flipflop.time guppy_basecaller --config dna_r9.4.1_450bps_flipflop.cfg --device auto -i $fast5_dir -t 12 -s $out_dir
    cat "$out_dir"/*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > guppy_v"$v"_flipflop.fastq
done

# For only the most recent version of Guppy, I also tried our custom-trained models:
v=2.2.3
out_dir="$output_prefix"/guppy_v"$v"_custom_kp_model
/usr/bin/time -v -o guppy_v"$v"_custom_kp_model.time guppy_basecaller --model_file holtlab_kp_r9.4_r9.4.1_nov_2018.jsn --config dna_r9.4.1_450bps.cfg --device auto -i $fast5_dir -t 12 -s $out_dir
cat "$out_dir"/*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > guppy_v"$v"_custom_kp_model.fastq
out_dir="$output_prefix"/guppy_v"$v"_custom_kp_big_net_model
/usr/bin/time -v -o guppy_v"$v"_custom_kp_big_net_model.time guppy_basecaller --model_file holtlab_kp_big_r9.4_r9.4.1_nov_2018.jsn --config dna_r9.4.1_450bps.cfg -i $fast5_dir -t 12 -s $out_dir
cat "$out_dir"/*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > guppy_v"$v"_custom_kp_big_net_model.fastq


for v in 1.5.1 1.6.0 1.8.1 1.8.3 1.8.5; do
    apt install "$guppy_installers"/ont-guppy_"$v"-1~xenial_amd64.deb
    out_dir="$output_prefix"/guppy_v"$v"
    /usr/bin/time -v -o guppy_v"$v".time guppy_basecaller --config dna_r9.4_450bps.cfg --device auto -i $fast5_dir -t 12 -s $out_dir
    cat "$out_dir"/*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > guppy_v"$v".fastq
done

# I couldn't make v1.4.3 run on the GPU so I used the CPU instead:
v=1.4.3
apt install "$guppy_installers"/ont-guppy_"$v"-1~xenial_amd64.deb
out_dir="$output_prefix"/guppy_v"$v"
/usr/bin/time -v -o guppy_v"$v".time guppy_basecaller --config dna_r9.4_450bps.cfg -i $fast5_dir -t 12 -s $out_dir
cat "$out_dir"/*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > guppy_v"$v".fastq

v=0.5.4
apt install "$guppy_installers"/ont-guppy_"$v"-1~xenial_amd64.deb
out_dir="$output_prefix"/guppy_v"$v"
/usr/bin/time -v -o guppy_v"$v".time guppy_basecaller --config dna_r9.4_450bps.cfg --device cuda:00000000:65:00.0 -i $fast5_dir -t 12 -s $out_dir
cat "$out_dir"/*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > guppy_v"$v".fastq

# I couldn't make v0.3.0 run on the GPU so I used the CPU instead:
v=0.3.0
apt install "$guppy_installers"/ont-guppy_"$v"-1~xenial_amd64.deb
out_dir="$output_prefix"/guppy_v"$v"
/usr/bin/time -v -o guppy_v"$v".time guppy --config dna_r9.4_450bps.cfg --device cpu -i $fast5_dir -s $out_dir
cat "$out_dir"/*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > guppy_v"$v".fastq
