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

# These commands assume you are running in a directory with a subdirectory titled
# 01_original_fast5s, which contains the isolate's fast5s. It will make a few more directories as
# it runs.


# Edit the following paths before running, as appropriate for your environment:
script_dir=/path/to/python_script_dir  # directory with other Python scripts (trim_signal.py, filter_reads.py and subdivide_read_dir.py)
r1=/path/to/reads_1.fastq.gz           # Illumina reads (pair 1)
r2=/path/to/reads_2.fastq.gz           # Illumina reads (pair 2)
sloika_dir=/path/to/sloika_dir         # directory of Sloika clone (https://github.com/rrwick/sloika)


# Make reference FASTA:
skesa --fastq "$r1","$r2" --cores 16 --memory 128 > ref_contigs.fasta

# Signal-level trimming of fast5s:
python3 "$script_dir"/trim_signal.py --trim_amount 2000 --min_size 50000 01_original_fast5s 02_trimmed_fast5s

# Basecalling and alignment - prep for QC:
read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -i 02_trimmed_fast5s -t 12 -s 03_basecalling -o fastq --disable_filtering --disable_pings
cat 03_basecalling/workspace/*.fastq > temp.fastq
minimap2 -c -x map-ont -t 20 ref_contigs.fasta temp.fastq > alignments.paf
rm temp.fastq

# Alignment-based QC and per-read references:
python3 "$script_dir"/filter_reads.py --min_basecalled_length 5000 --max_unaligned_bases 30 --max_window_indels 0.8 --window_size 25 02_trimmed_fast5s 03_basecalling/sequencing_summary.txt ref_contigs.fasta alignments.paf 04_filtered_fast5s read_references.fasta

# Run the first round of Sloika chunky in parallel on smallish groups of reads. This means that
# troublesome reads which hang and/or use lots of RAM only mess up their batch, not the whole
# group.
python3 "$script_dir"/subdivide_read_dir.py 04_filtered_fast5s 50
cd 04_filtered_fast5s
for d in *; do
    if (( d < 240 )); then  # limit to 240 batches (12k reads)
        sbatch --job-name=Sloika_chunkify_"$d" --account=js66 --ntasks=4 --nodes=1-1 --time=1:00:00 --mem=40000 --wrap "module load python/3.5.2-gcc5; source ~/programs/sloika/build/env/bin/activate; export THEANO_FLAGS=device=cpu,floatX=float32,mode=FAST_RUN,blas.ldflags='-lblas',scan.allow_gc=False; /usr/bin/time -v ~/programs/sloika/bin/chunkify.py raw_remap --overwrite --jobs 4 --chunk_len 4000 --downsample_factor 5 --output_strand_list ../unfiltered_strands.txt "$d" remapped_unfiltered_"$d".hdf5 ~/programs/sloika/models/pretrained.pkl ../read_references.fasta"
        sleep 0.1
    fi
done
rm remapped_unfiltered_*.hdf5
cd ..

# Optionally, visualise the distributions of scores and proportions of stays:
tail -n+2 unfiltered_strands.txt | awk '{print $3;}' | histogram.py --min 0 --max 1.6 --buckets 40
tail -n+2 unfiltered_strands.txt | awk '{print $4 / ($7 - $6 + $4);}' | histogram.py --min 0 --max 1 --buckets 40

# QC pass using the chunkify strand info:
low_score=$(tail -n+2 unfiltered_strands.txt | awk '{print $3;}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.15 - 0.5)]}')
high_score=$(tail -n+2 unfiltered_strands.txt | awk '{print $3;}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
min_cov=$(tail -n+2 unfiltered_strands.txt | awk '{print ($7 - $6) / $5;}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.1 - 0.5)]}')
min_cov=$(( min_cov > 0.975 ? min_cov : 0.975 ))
max_stay=$(tail -n+2 unfiltered_strands.txt | awk '{print $4 / ($7 - $6 + $4);}' | sort -n | awk '{all[NR] = $0} END{print all[int(NR*0.75 - 0.5)]}')
(head -n 1 unfiltered_strands.txt ; cat unfiltered_strands.txt | awk -v a="$low_score" -v b="$high_score" -v c="$min_cov" -v d="$max_stay" '$3 > a && $3 < b && ($7 - $6) > c * $5 && $4 / ($7 - $6 + $4) < d && $6 > 0') > filtered_strand_list.txt

# Split reads into a training set (most of them) and a validation set (a small fraction):
mkdir 05_training_fast5s
for f in $(tail -n+2 filtered_strand_list.txt | cut -f1); do
    cp 04_filtered_fast5s/*/"$f" 05_training_fast5s
done
mkdir 06_validation_fast5s
for f in $(ls 05_training_fast5s | shuf | head -n 20); do
    mv 05_training_fast5s/"$f" 06_validation_fast5s
done
python3 "$script_dir"/subdivide_read_dir.py 05_training_fast5s 1000

# The second round of Sloika chunkify:
isolate_name=$(basename "$PWD")
cd 05_training_fast5s
for d in *; do
    sbatch --job-name=Chunkify_"$isolate_name"_"$d" --account=js66 --ntasks=12 --nodes=1-1 --time=24:00:00 --mem=120000 --wrap "module load python/3.5.2-gcc5; source ~/programs/sloika/build/env/bin/activate; export THEANO_FLAGS=device=cpu,floatX=float32,mode=FAST_RUN,blas.ldflags='-lblas',scan.allow_gc=False; /usr/bin/time -v ~/programs/sloika/bin/chunkify.py raw_remap --jobs 12 --chunk_len 4000 --downsample_factor 5 --output_strand_list strands_"$d".txt "$d" remapped_"$d".hdf5 ~/programs/sloika/models/pretrained.pkl ../read_references.fasta"
    sleep 0.1
done
