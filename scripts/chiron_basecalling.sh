#!/usr/bin/env bash

fast5_dir=/var/lib/MinKNOW/data/basecalling_comparison/fast5s
output_prefix=/var/lib/MinKNOW/data/basecalling_comparison/basecalling_output
chiron_dirs=/home/holt-ont/Chiron

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
        /usr/bin/time -v -o $time_file python "$chiron_dirs"/Chiron-"$v"/chiron/entry.py call -i ../fast5s_100_read_dirs/"$d" -o $out_dir || break
        cat "$out_dir"/result/*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > $out_file
    fi
done


# Older versions (v0.2 and v0.3) need TensorFlow 1.0.1, which needs CUDA 8.0.
source /home/holt-ont/Chiron/chiron_virtualenv/bin/activate
export PATH=/home/holt-ont/Chiron/cuda-8.0/bin:$PATH
export LD_LIBRARY_PATH=/home/holt-ont/Chiron/cuda-8.0/lib64:$LD_LIBRARY_PATH

for v in 0.3 0.2; do
    for d in {000..151}; do
        out_dir="$output_prefix"/chiron_v"$v"_"$d"
        out_file=chiron_v"$v"_"$d".fastq
        time_file=chiron_v"$v"_"$d".time
        if [ ! -f $out_file ]; then
            rm -rf $out_dir
            rm -f $time_file
            /usr/bin/time -v -o $time_file python "$chiron_dirs"/Chiron-"$v"/chiron/entry.py call -i ../fast5s_100_read_dirs/"$d" -o $out_dir || break
            cat "$out_dir"/result/*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > $out_file
        fi
    done
done

cat chiron_v0.2_*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > chiron_v0.2.fastq
cat chiron_v0.3_*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > chiron_v0.3.fastq
cat chiron_v0.4.2_*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > chiron_v0.4.2.fastq

cat chiron_v0.2_*.time > chiron_v0.2.time
cat chiron_v0.3_*.time > chiron_v0.3.time
cat chiron_v0.4.2_*.time > chiron_v0.4.2.time
