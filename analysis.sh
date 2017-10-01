#!/usr/bin/env bash

# This script conducts read and assembly analysis on a set of ONT reads, comparing them to a
# reference sequence. It is part of https://github.com/rrwick/Basecalling-comparison.

# This script expects to find the following in the directory where it's run:
#   * reference.fasta: the reference sequence
#   * 01_raw_fast5 directory: has all fast5 files
#   * 02_basecalled_reads directory: has one or more fastq.gz/fasta.gz read files
#   * read_id_to_fast5: a file with two columns: read_ID in the first and fast5 filename in the second
#   * illumina_1.fastq.gz illumina_2.fastq.gz: Illumina reads for the same sample

# Set this to the desired number of threads (for alignment and polishing).
threads=32

# Set this to the directory containing the Python scripts (e.g. read_length_identity.py).
python_script_dir=/path/to/python_scripts

# Set the full path to Nanopolish here.
nanopolish_exec_dir=/path/to/nanopolish

# If you want to run this script on all read files in the 02_basecalled_reads directory, leave
# these lines uncommented:
cd 02_basecalled_reads
read_files=$(ls)
cd ..

# If you want to run this script only on particular read files, change and uncomment this line:
# read_files="albacore_v0.8.4.fastq.gz"

# Make necessary directories.
mkdir -p 03_read_names_fixed
mkdir -p 04_read_data
mkdir -p 05_trimmed_reads
mkdir -p 06_subsampled_reads
mkdir -p 07_assemblies
mkdir -p 08_assembly_data
mkdir -p 09_nanopolish
mkdir -p 10_nanopolish_data

# Create a table of basic info about each read.
python3 "$python_script_dir"/read_table.py 01_raw_fast5 > 04_read_data/read_data.tsv

for f in $read_files; do
    set=${f%.fastq.gz}
    set=${set%.fasta.gz}

    # Save file paths in variables, for brevity.
    raw_fast5_dir=01_raw_fast5
    all_reads=02_basecalled_reads/"$f"
    all_reads_fixed_names=03_read_names_fixed/"$f"
    read_alignment=04_read_data/"$set".paf
    read_data=04_read_data/"$set"_reads.tsv
    trimmed_reads=05_trimmed_reads/"$set".fastq.gz
    subsampled_reads=06_subsampled_reads/"$set".fastq.gz
    assembly=07_assemblies/"$set"_assembly.fasta
    assembly_pieces=08_assembly_data/"$set"_assembly_pieces.fasta
    assembly_alignment=08_assembly_data/"$set"_assembly.paf
    assembly_data=08_assembly_data/"$set"_assembly.tsv
    nanopolish_assembly_dir=09_nanopolish
    nanopolish_assembly=09_nanopolish/"$set"_nanopolish.fasta
    nanopolish_assembly_pieces=10_nanopolish_data/"$set"_nanopolish_pieces.fasta
    nanopolish_assembly_alignment=10_nanopolish_data/"$set"_nanopolish.paf
    nanopolish_assembly_data=10_nanopolish_data/"$set"_nanopolish.tsv

    printf "\n\n\n\n"
    echo "NORMALISE READ HEADERS: "$set
    echo "--------------------------------------------------------------------------------"
    python3 "$python_script_dir"/fix_read_names.py $all_reads read_id_to_fast5 | gzip > $all_reads_fixed_names

    printf "\n\n\n\n"
    echo "ASSESS READS: "$set
    echo "--------------------------------------------------------------------------------"
    minimap2 -k12 -t $threads -c reference.fasta $all_reads_fixed_names > $read_alignment
    python3 "$python_script_dir"/read_length_identity.py $all_reads_fixed_names $read_alignment > $read_data
    rm $read_alignment

    printf "\n\n\n\n"
    echo "ASSEMBLY: "$set
    echo "--------------------------------------------------------------------------------"
    porechop -i $all_reads_fixed_names -o $trimmed_reads --no_split --threads $threads --check_reads 1000
    filtlong -1 illumina_1.fastq.gz -2 illumina_2.fastq.gz --min_length 1000 --target_bases 500000000 --trim --split 250 $trimmed_reads | gzip > $subsampled_reads
    rebaler -t $threads reference.fasta $subsampled_reads > $assembly

    printf "\n\n\n\n"
    echo "ASSESS ASSEMBLY: "$set
    echo "--------------------------------------------------------------------------------"
    python3 "$python_script_dir"/chop_up_assembly.py $assembly 10000 > $assembly_pieces
    minimap2 -k12 -t $threads -c reference.fasta $assembly_pieces > $assembly_alignment
    python3 ../read_length_identity.py $assembly_pieces $assembly_alignment > $assembly_data
    rm $assembly_pieces $assembly_alignment

    printf "\n\n\n\n"
    echo "NANOPOLISH: "$set
    echo "--------------------------------------------------------------------------------"
    python3 "$python_script_dir"/nanopolish_slurm_wrapper.py $assembly $all_reads_fixed_names $raw_fast5_dir $nanopolish_assembly_dir $nanopolish_exec_dir $threads
    rm "$all_reads_fixed_names".fa.gz*
    rm "$assembly".fai

    printf "\n\n\n\n"
    echo "ASSESS NANOPOLISHED ASSEMBLY: "$set
    echo "--------------------------------------------------------------------------------"
    python3 "$python_script_dir"/chop_up_assembly.py $nanopolish_assembly 10000 > $nanopolish_assembly_pieces
    minimap2 -x map10k -t $threads -c reference.fasta $nanopolish_assembly_pieces > $nanopolish_assembly_alignment
    python3 "$python_script_dir"/read_length_identity.py $nanopolish_assembly_pieces $nanopolish_assembly_alignment > $nanopolish_assembly_data
    rm $nanopolish_assembly_pieces $nanopolish_assembly_alignment

done
