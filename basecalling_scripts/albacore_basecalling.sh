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

# This script will execute each version of Albacore on the reads.


# Edit the following paths before running, as appropriate for your environment:
fast5_dir=/path/to/fast5s                         # all fast5s must be directly in this directory (i.e. not nested in subdirectories)
output_prefix=/path/to/basecalling_output         # Albacore output directories and consolidated fastqs will go in this directory
albacore_installers=/path/to/albacore_installers  # must contain the whl/deb files for the Albacore versions - this script installs each version before running it
fast5_to_fastq=/path/to/fastq_to_fastq.py         # path to fastq_to_fastq.py script (from https://github.com/rrwick/Fast5-to-Fastq) - required for older versions of Albacore


cd $output_prefix

# The most recent versions of Albacore:
for v in 2.1.0 2.1.1 2.1.1 2.1.2 2.1.3 2.1.5 2.1.6 2.1.7 2.1.9 2.1.10 2.2.1 2.2.2 2.2.4 2.2.6 2.2.7 2.3.0 2.3.1 2.3.2 2.3.3 2.3.4; do
    pip3 install "$albacore_installers"/ont_albacore-"$v"-cp35-cp35m-manylinux1_x86_64.whl
    out_dir="$output_prefix"/albacore_v"$v"
    /usr/bin/time -v -o albacore_v"$v".time read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -i $fast5_dir -t 12 -s $out_dir -o fastq --disable_filtering --disable_pings
    cat "$out_dir"/workspace/*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > albacore_v"$v".fastq
done

# These versions don't have --disable_pings:
for v in 2.0.0 2.0.1 2.0.2; do
    pip3 install "$albacore_installers"/ont_albacore-"$v"-cp35-cp35m-manylinux1_x86_64.whl
    out_dir="$output_prefix"/albacore_v"$v"
    /usr/bin/time -v -o albacore_v"$v".time read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -i $fast5_dir -t 12 -s $out_dir -o fastq --disable_filtering
    cat "$out_dir"/workspace/*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > albacore_v"$v".fastq
done

# These versions don't have --disable_filtering:
for v in 1.2.0 1.2.1 1.2.2 1.2.3 1.2.4 1.2.5 1.2.6 1.3.1; do
    pip3 install "$albacore_installers"/ont_albacore-"$v"-cp35-cp35m-manylinux1_x86_64.whl
    out_dir="$output_prefix"/albacore_v"$v"
    /usr/bin/time -v -o albacore_v"$v".time read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -i $fast5_dir -t 12 -s $out_dir -o fastq
    cat "$out_dir"/workspace/*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > albacore_v"$v".fastq
done

# These versions mess up the fastq files (extra line break) which requires fixing:
for v in 1.1.0 1.1.1 1.1.2; do
    pip3 install "$albacore_installers"/ont_albacore-"$v"-cp35-cp35m-manylinux1_x86_64.whl
    out_dir="$output_prefix"/albacore_v"$v"
    /usr/bin/time -v -o albacore_v"$v".time read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -i $fast5_dir -t 12 -s $out_dir -o fastq
    cat "$out_dir"/workspace/*.fastq | paste - - - - - | sort -k1,1 -t " " | tr "\t" "\n" > albacore_v"$v".fastq
done

# These versions can't output to fastq:
for v in 1.0.1 1.0.2 1.0.3 1.0.4; do
    pip3 install "$albacore_installers"/ont_albacore-"$v"-cp35-cp35m-manylinux1_x86_64.whl
    out_dir="$output_prefix"/albacore_v"$v"
    /usr/bin/time -v -o albacore_v"$v".time read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -i $fast5_dir -t 12 -s $out_dir
    ~/Fast5-to-Fastq/fast5_to_fastq.py "$out_dir"/workspace | paste - - - - | sed 's/_Basecall_1D_template/ Basecall_1D_template/' | sort -k1,1 -t " " | tr "\t" "\n" > albacore_v"$v".fastq
done

# These versions don't have flowcell and kit options (need to specify config):
for v in 0.8.2 0.8.4 0.9.1; do
    pip3 install "$albacore_installers"/ont_albacore-"$v"-cp35-cp35m-manylinux1_x86_64.whl
    out_dir="$output_prefix"/albacore_v"$v"
    /usr/bin/time -v -o albacore_v"$v".time read_fast5_basecaller.py -c FLO-MIN106_LSK108_linear.cfg -i $fast5_dir -t 12 -s $out_dir
    "$fast5_to_fastq" "$out_dir"/workspace | paste - - - - | sed 's/_Basecall_1D_template/ Basecall_1D_template/' | sort -k1,1 -t " " | tr "\t" "\n" > albacore_v"$v".fastq
done

# 0.7.5 is the only version without a Python wheel:
v=0.7.5
apt install "$albacore_installers"/python3-ont-albacore_"$v"-1-xenial_amd64.deb
out_dir="$output_prefix"/albacore_v"$v"
/usr/bin/time -v -o albacore_v"$v".time read_fast5_basecaller.py -c FLO-MIN106_LSK108_linear.cfg -i $fast5_dir -t 12 -s $out_dir
"$fast5_to_fastq" "$out_dir"/workspace | paste - - - - | sed 's/_Basecall_1D_template/ Basecall_1D_template/' | sort -k1,1 -t " " | tr "\t" "\n" > albacore_v"$v".fastq
