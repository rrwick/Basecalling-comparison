#!/bin/bash
set -e

# Starting with the binned set of reads for this isolate, copy all of the corresponding fast5 files to a local directory.
cd /home/UNIMELB/inouye-hpc-sa/nanopore-data/albacore_comparison
zcat /MDHS/Research/SysGen-Lab/MinION/2017-04-24_Kleb_barcode/nanopore_reads/fastq/2_trimmed/barcode08/BC08.fastq.gz | grep "^@" | grep -P -o "\w{8}-\w{4}-\w{4}-\w{4}-\w{12}" > read_ids
grep -F -f read_ids /MDHS/Research/SysGen-Lab/MinION/2017-04-24_Kleb_barcode/nanopore_reads/basecalling/to_fastq/sequencing_summary.txt | grep -o -P "[\w_]+\.fast5" > read_filenames
find /MDHS/Research/SysGen-Lab/MinION/2017-04-24_Kleb_barcode/nanopore_reads/raw_fast5 -name "*.fast5" > all_read_filenames
grep -F -f read_filenames all_read_filenames > read_filenames_with_path
mkdir raw_fast5
for read in $(cat read_filenames_with_path); do cp $read raw_fast5; done
mkdir gather_reads
albacore_comparison mv all_read_filenames read_filenames read_filenames_with_path read_ids read_names gather_reads


assembly_identity_distribution() {
    python3 chop_up_assembly.py "$1" 10000 > "$1".pieces.fasta
    minimap2 -x map10k -t $threads -c Kleb_KSB1_7J.fasta "$1".fastq.gz > "$1".paf
    python3 read_length_identity.py "$1".paf > "$1".tsv
    python3 histograms.py "$1".tsv 0.1 0.1 "$1".identity_histogram "$1".relative_length_histogram
}


extract_map_and_assemble () {
    threads=32

    # Extract fastq from fast5 files
    ~/Fast5-to-Fastq/fast5_to_fastq.py "$1" | gzip > "$1".fastq.gz

    # Align reads to reference and get reads identities
    minimap2 -x map10k -t $threads -c Kleb_KSB1_7J.fasta "$1".fastq.gz > albacore_v2.0.0.paf
    python3 read_length_identity.py "$1".paf > "$1"_reads.tsv
    (head -n 1 "$1".tsv && tail -n +2 "$1"_reads.tsv | sort) > "$1"_reads_sorted.tsv; mv "$1"_reads_sorted.tsv "$1"_reads.tsv

    # Trim and subsample reads for assembly
    porechop -i "$1".fastq.gz -o "$1"_trimmed.fastq.gz --no_split --threads $threads --check_reads 1000
    filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 "$1"_trimmed.fastq.gz | gzip > "$1"_subsampled.fastq.gz

    # Do a Nanopore-only assembly and get assembly identities
    unicycler -l "$1"_subsampled.fastq.gz -o "$1"_assembly --threads $threads
    assembly_identity_distribution("$1"_assembly/assembly.fasta)

    # Improve the assembly with Nanopolish and get assembly identities again
    mkdir "$1"_nanopolish
    nanopolish extract --type template "$1" > "$1"_nanopolish/reads.fa
    bwa mem -x ont2d -t $threads "$1"_assembly/assembly.fasta "$1"_nanopolish/reads.fa | samtools sort -o "$1"_nanopolish/reads.sorted.bam -T reads.tmp -
    samtools index "$1"_nanopolish/reads.sorted.bam
    python nanopolish_makerange.py "$1"_assembly/assembly.fasta | parallel --results nanopolish.results -P 8 nanopolish variants --consensus "$1"_nanopolish/polished.{1}.fa -w {1} -r "$1"_nanopolish/reads.fa -b "$1"_nanopolish/reads.sorted.bam -g "$1"_assembly/assembly.fasta -t 4 --min-candidate-frequency 0.1
    python nanopolish_merge.py "$1"_nanopolish/polished.*.fa > "$1"_nanopolish/polished_genome.fa
    assembly_identity_distribution("$1"_nanopolish/assembly.fasta)
}


# Albacore v2.0.0
cd /home/UNIMELB/inouye-hpc-sa
pip3 install ont_albacore-2.0.0-cp35-cp35m-manylinux1_x86_64.whl
cd /home/UNIMELB/inouye-hpc-sa/nanopore-data/albacore_comparison
read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -i raw_fast5 -t $threads -s albacore_v2.0.0 -o fast5 --disable_filtering
extract_map_and_assemble("albacore_v2.0.0")


# Albacore v1.2.6
cd /home/UNIMELB/inouye-hpc-sa
pip3 install ont_albacore-1.2.6-cp35-cp35m-manylinux1_x86_64.whl
cd /home/UNIMELB/inouye-hpc-sa/nanopore-data/albacore_comparison
read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -i raw_fast5 -t $threads -s albacore_v1.2.6 -o fast5
extract_map_and_assemble("albacore_v1.2.6")


# Albacore v1.1.2
cd /home/UNIMELB/inouye-hpc-sa
pip3 install ont_albacore-1.1.2-cp35-cp35m-manylinux1_x86_64.whl
cd /home/UNIMELB/inouye-hpc-sa/nanopore-data/albacore_comparison
read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -i raw_fast5 -t $threads -s albacore_v1.1.2 -o fast5
extract_map_and_assemble("albacore_v1.1.2")


# Albacore v1.0.4
cd /home/UNIMELB/inouye-hpc-sa
pip3 install ont_albacore-1.0.4-cp35-cp35m-manylinux1_x86_64.whl
cd /home/UNIMELB/inouye-hpc-sa/nanopore-data/albacore_comparison
read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -i raw_fast5 -t $threads -s albacore_v1.0.4
extract_map_and_assemble("albacore_v1.0.4")


# Albacore v0.9.1
cd /home/UNIMELB/inouye-hpc-sa
pip3 install ont_albacore-0.9.1-cp35-cp35m-manylinux1_x86_64.whl
cd /home/UNIMELB/inouye-hpc-sa/nanopore-data/albacore_comparison
read_fast5_basecaller.py -c FLO-MIN106_LSK108_linear.cfg -i raw_fast5 -t $threads -s albacore_v0.9.1
extract_map_and_assemble("albacore_v0.9.1")


# Albacore v0.8.4
cd /home/UNIMELB/inouye-hpc-sa
pip3 install ont_albacore-0.8.4-cp35-cp35m-manylinux1_x86_64.whl
cd /home/UNIMELB/inouye-hpc-sa/nanopore-data/albacore_comparison
read_fast5_basecaller.py -c FLO-MIN106_LSK108_linear.cfg -i raw_fast5 -t $threads -s albacore_v0.8.4
extract_map_and_assemble("albacore_v0.8.4")
