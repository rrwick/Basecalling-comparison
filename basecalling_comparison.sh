#!/bin/bash


gather_fast5s_files=false
albacore_v2_0_1=false
albacore_v2_0_0=false
albacore_v1_2_6=false
albacore_v1_1_2=false
albacore_v1_0_4=false
albacore_v0_9_1=false
albacore_v0_8_4=false
scrappie_v1_1_0_raw=false
scrappie_v1_1_0_events=false
scrappie_v1_0_0_raw=false
scrappie_v1_0_0_events=false
nanonet=false
chiron=false


threads=40


# Local paths specific to this server.
raw_fast5=/MDHS/Research/SysGen-Lab/MinION/2017-04-24_Kleb_barcode/nanopore_reads/raw_fast5
all_reads=/MDHS/Research/SysGen-Lab/MinION/2017-04-24_Kleb_barcode/nanopore_reads/fastq/2_trimmed/barcode01/BC01.fastq.gz
albacore_table=/MDHS/Research/SysGen-Lab/MinION/2017-04-24_Kleb_barcode/nanopore_reads/basecalling/to_fastq/sequencing_summary.txt
reference=/MDHS/Research/SysGen-Lab/Long_read_assembly/Kleb_INF042/Kleb_INF042.fasta
albacore_whl_dir=/home/UNIMELB/inouye-hpc-sa
nanopolish_scripts_dir=/home/UNIMELB/inouye-hpc-sa/nanopolish/scripts
scrappie_v1_0_0_path=/home/UNIMELB/inouye-hpc-sa/scrappie-release-1.0.0/build
scrappie_v1_1_0_path=/home/UNIMELB/inouye-hpc-sa/scrappie-release-1.1.0/build


if $gather_fast5s_files; then

    # Starting with the binned set of reads for this isolate, copy all of the corresponding fast5
    # files to a local directory and then delete the smaller ones.
    mkdir gather_reads
    zgrep "^@" $all_reads | grep -P -o "\w{8}-\w{4}-\w{4}-\w{4}-\w{12}" > gather_reads/read_ids
    grep -F -f gather_reads/read_ids $albacore_table | grep -o -P "[\w_]+\.fast5" > gather_reads/read_filenames
    find $raw_fast5 -name "*.fast5" > gather_reads/all_read_filenames
    grep -F -f gather_reads/read_filenames gather_reads/all_read_filenames > gather_reads/read_filenames_with_path
    mkdir raw_fast5
    for read in $(cat gather_reads/read_filenames_with_path); do cp $read raw_fast5; done
    find raw_fast5 -name "*.fast5" -size -100k -delete
    rm -r gather_reads
fi


# This is the reference used to assess the reads and assemblies.
cp $reference reference.fasta
cp $nanopolish_scripts_dir/nanopolish_makerange.py .
cp $nanopolish_scripts_dir/nanopolish_merge.py .


assembly_identity_distribution() {
    python3 chop_up_assembly.py "$1" 10000 > "$2"_pieces.fasta
    minimap2 -x map10k -t $threads -c reference.fasta "$2"_pieces.fasta > "$2".paf
    python3 read_length_identity.py "$2".paf > "$2".tsv
    rm "$2"_pieces.fasta
}


extract_map_and_assemble () {

     # Some basecallers (e.g. Albacore) produce a directory of fast5s from which we need to extract
     # a fastq. Others (e.g. Chiron and Scrappie) just produce a fasta.
    if [ -d "$1" ]; then
        fast5_to_fastq.py "$1"/workspace | gzip > "$1".fastq.gz
        basecalled_reads="$1".fastq.gz
    else
        basecalled_reads="$1".fasta
    fi

    # Align reads to reference and get reads identities
    minimap2 -x map10k -t $threads -c reference.fasta $basecalled_reads > "$1".paf
    python3 read_length_identity.py "$1".paf > "$1"_reads.tsv
    (head -n 1 "$1"_reads.tsv && tail -n +2 "$1"_reads.tsv | sort) > "$1"_reads_sorted.tsv; mv "$1"_reads_sorted.tsv "$1"_reads.tsv
    python3 histograms.py "$1"_reads.tsv 0.1 0.1 "$1"_read_identity_histogram "$1"_read_relative_length_histogram

    # Trim and subsample reads for assembly
    porechop -i $basecalled_reads -o "$1"_trimmed.fastq.gz --no_split --threads $threads --check_reads 1000
    filtlong -1 illumina_1.fastq.gz -2 illumina_2.fastq.gz --min_length 1000 --target_bases 500000000 "$1"_trimmed.fastq.gz | gzip > "$1"_subsampled.fastq.gz

    # Do a Nanopore-only assembly and get assembly identities
    unicycler -l "$1"_subsampled.fastq.gz -o "$1"_assembly --threads $threads
    cp "$1"_assembly/assembly.fasta "$1"_assembly.fasta
    assembly_identity_distribution "$1"_assembly.fasta "$1"_assembly

    # Improve the assembly with Nanopolish and get assembly identities again
    mkdir "$1"_nanopolish
    nanopolish extract --type template --recurse "$1"/workspace > "$1"_nanopolish/reads.fa
    bwa index "$1"_assembly.fasta
    bwa mem -x ont2d -t $threads "$1"_assembly.fasta "$1"_nanopolish/reads.fa | samtools sort -o "$1"_nanopolish/reads.sorted.bam -T reads.tmp -
    rm "$1"_assembly.fasta.*  # clean up BWA index files
    samtools index "$1"_nanopolish/reads.sorted.bam
    python nanopolish_makerange.py "$1"_assembly.fasta | parallel --results "$1"_nanopolish/nanopolish.results -P 10 nanopolish variants --consensus "$1"_nanopolish/polished.{1}.fa -w {1} -r "$1"_nanopolish/reads.fa -b "$1"_nanopolish/reads.sorted.bam -g "$1"_assembly.fasta -t 4 --min-candidate-frequency 0.1
    python nanopolish_merge.py "$1"_nanopolish/polished.*.fa > "$1"_nanopolished_assembly.fasta
    assembly_identity_distribution "$1"_nanopolished_assembly.fasta "$1"_nanopolished_assembly
}


if $albacore_v2_0_1; then
    pip3 install $albacore_whl_dir/ont_albacore-2.0.1-cp35-cp35m-manylinux1_x86_64.whl
    read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -i raw_fast5 -t $threads -s albacore_v2.0.1 -o fast5 --disable_filtering
    extract_map_and_assemble "albacore_v2.0.1"
fi

if $albacore_v2_0_0; then
    pip3 install $albacore_whl_dir/ont_albacore-2.0.0-cp35-cp35m-manylinux1_x86_64.whl
    read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -i raw_fast5 -t $threads -s albacore_v2.0.0 -o fast5 --disable_filtering
    extract_map_and_assemble "albacore_v2.0.0"
fi

if $albacore_v1_2_6; then
    pip3 install $albacore_whl_dir/ont_albacore-1.2.6-cp35-cp35m-manylinux1_x86_64.whl
    read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -i raw_fast5 -t $threads -s albacore_v1.2.6 -o fast5
    extract_map_and_assemble "albacore_v1.2.6"
fi

if $albacore_v1_1_2; then
    pip3 install $albacore_whl_dir/ont_albacore-1.1.2-cp35-cp35m-manylinux1_x86_64.whl
    read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -i raw_fast5 -t $threads -s albacore_v1.1.2 -o fast5
    extract_map_and_assemble "albacore_v1.1.2"
fi

if $albacore_v1_0_4; then
    pip3 install $albacore_whl_dir/ont_albacore-1.0.4-cp35-cp35m-manylinux1_x86_64.whl
    read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -i raw_fast5 -t $threads -s albacore_v1.0.4
    extract_map_and_assemble "albacore_v1.0.4"
fi

if $albacore_v0_9_1; then
    pip3 install $albacore_whl_dir/ont_albacore-0.9.1-cp35-cp35m-manylinux1_x86_64.whl
    read_fast5_basecaller.py -c FLO-MIN106_LSK108_linear.cfg -i raw_fast5 -t $threads -s albacore_v0.9.1
    extract_map_and_assemble "albacore_v0.9.1"
fi

if $albacore_v0_8_4; then
    pip3 install $albacore_whl_dir/ont_albacore-0.8.4-cp35-cp35m-manylinux1_x86_64.whl
    read_fast5_basecaller.py -c FLO-MIN106_LSK108_linear.cfg -i raw_fast5 -t $threads -s albacore_v0.8.4
    extract_map_and_assemble "albacore_v0.8.4"
fi

if $scrappie_v1_0_0_raw; then
    export OMP_NUM_THREADS=$threads
    export OPENBLAS_NUM_THREADS=1
    $scrappie_v1_0_0_path/scrappie raw raw_fast5 --threads=$threads > scrappie_v1_0_0_raw.fasta
    extract_map_and_assemble "scrappie_v1_0_0_raw"
fi

if $scrappie_v1_0_0_events; then
    export OMP_NUM_THREADS=$threads
    export OPENBLAS_NUM_THREADS=1
    mkdir albacore_v1.2.6_fast5s
    for f in $(find albacore_v1.2.6/workspace -name "*.fast5"); do cp $f albacore_v1.2.6_fast5s; done
    $scrappie_v1_0_0_path/scrappie events albacore_v1.2.6_fast5s --threads=$threads --albacore > scrappie_v1_0_0_events.fasta
    extract_map_and_assemble "scrappie_v1_0_0_events"
fi

if $scrappie_v1_1_0_raw; then
    export OMP_NUM_THREADS=$threads
    export OPENBLAS_NUM_THREADS=1
    $scrappie_v1_1_0_path/scrappie raw raw_fast5 --threads=$threads > scrappie_v1_1_0_raw.fasta
    extract_map_and_assemble "scrappie_v1_1_0_raw"
fi

if $scrappie_v1_1_0_events; then
    export OMP_NUM_THREADS=$threads
    export OPENBLAS_NUM_THREADS=1
    $scrappie_v1_1_0_path/scrappie events raw_fast5 --threads=$threads --albacore > scrappie_v1_1_0_events.fasta
    extract_map_and_assemble "scrappie_v1_1_0_events"
fi

if $nanonet; then
    cp -r raw_fast5 nanonet  # Nanonet adds data to the fast5s, so we first make a copy.
    nanonetcall --chemistry r9.4 --write_events --min_len 1 --max_len 1000000 --jobs $threads nanonet > /dev/null
    extract_map_and_assemble "nanonet"
fi

if $chiron; then
    source /home/UNIMELB/inouye-hpc-sa/chiron/chiron/bin/activate
    chiron call -i raw_fast5 -o chiron -t $threads
    paste --delimiter=\\n --serial chiron/result/*.fasta > chiron.fasta
    deactivate
    extract_map_and_assemble "chiron"
fi
