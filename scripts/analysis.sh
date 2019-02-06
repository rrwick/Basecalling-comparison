#!/usr/bin/env bash

# Copyright 2018 Ryan Wick (rrwick@gmail.com)
# https://github.com/rrwick/Basecalling-comparison

# This program is free software: you can redistribute it and/or modify it under the terms of the
# GNU General Public License as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version. This program is distributed in the hope that it
# will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You
# should have received a copy of the GNU General Public License along with this program. If not,
# see <http://www.gnu.org/licenses/>.

# This script:
#   * takes one argument: a fastq/fasta file of reads for analysis:
#      * read names should be in UUID format
#      * reads should be sorted by name
#   * expects the following in the directory in which it's run:
#      * a reference.fasta file
#   * expects the following tools in the PATH:
#      * minimap2
#      * rebaler
#      * racon
#      * nanopolish
#      * nanopolish_makerange.py
#      * MUMmer tools (nucmer, delta-filter, show-snps)


# Some high-level settings for the script:
threads=20
assembly_dir=assemblies     # where finished assemblies will go
nanopolish_dir=nanopolish   # where Nanopolished assemblies will go
results_dir=results         # where the result tables will go


# This variable holds the directory of the other scripts (assumed to be in the same directory as
# this script):
scripts_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"


# Get the read set name from the filename:
set=$(basename $1 | sed 's/.fastq.gz//' | sed 's/.fasta.gz//' | sed 's/.fastq//' | sed 's/.fasta//')
if [[ $1 = *"fastq"* ]]; then read_type=fastq
elif [[ $1 = *"fasta"* ]]; then read_type=fasta
else echo "Error: cannot determine read type"; exit 1
fi


# Various file paths that will be used along the way:
read_data="$results_dir"/"$set"_reads.tsv
read_alignment="$results_dir"/"$set"_reads.paf
trimmed_reads="$assembly_dir"/"$set"_trimmed."$read_type"
filtlong_reads="$assembly_dir"/"$set"_filtlong."$read_type"
assembly_reads="$assembly_dir"/"$set"_assembly_reads.fastq
final_assembly="$assembly_dir"/"$set"_assembly.fasta
assembly_pieces="$results_dir"/"$set"_assembly_pieces.fasta
assembly_data="$results_dir"/"$set"_assembly.tsv
assembly_alignment="$results_dir"/"$set"_assembly.paf
nanopolish_assembly="$nanopolish_dir"/"$set"_nanopolish.fasta
sequencing_summary=03_basecalling_output/"$set"/sequencing_summary.txt
bam_file="$1".sorted.bam
nanopolish_pieces="$results_dir"/"$set"_nanopolish_pieces.fasta
nanopolish_data="$results_dir"/"$set"_nanopolish.tsv
nanopolish_alignment="$results_dir"/"$set"_nanopolish.paf
nanopolish_2_assembly="$nanopolish_dir"/"$set"_nanopolish_2.fasta
bam_file_2="$1"_2.sorted.bam
nanopolish_2_pieces="$results_dir"/"$set"_nanopolish_2_pieces.fasta
nanopolish_2_data="$results_dir"/"$set"_nanopolish_2.tsv
nanopolish_2_alignment="$results_dir"/"$set"_nanopolish_2.paf
nanopolish_3_assembly="$nanopolish_dir"/"$set"_nanopolish_3.fasta
bam_file_3="$1"_3.sorted.bam
nanopolish_3_pieces="$results_dir"/"$set"_nanopolish_3_pieces.fasta
nanopolish_3_data="$results_dir"/"$set"_nanopolish_3.tsv
nanopolish_3_alignment="$results_dir"/"$set"_nanopolish_3.paf
nanopolish_4_assembly="$nanopolish_dir"/"$set"_nanopolish_4.fasta
bam_file_4="$1"_4.sorted.bam
nanopolish_4_pieces="$results_dir"/"$set"_nanopolish_4_pieces.fasta
nanopolish_4_data="$results_dir"/"$set"_nanopolish_4.tsv
nanopolish_4_alignment="$results_dir"/"$set"_nanopolish_4.paf


# Used for nucmer:
prefix="$set"_details
ref_contig=chromosome
assembly_contig=chromosome




printf "\n\n\n"
echo "ASSESS READS: "$set
echo "--------------------------------------------------------------------------------"
mkdir -p "$results_dir"
minimap2 -x map-ont -t $threads -c reference.fasta $1 > $read_alignment
pypy3 "$scripts_dir"/read_length_identity.py $1 $read_alignment > $read_data
rm $read_alignment




printf "\n\n\n\n"
echo "ASSEMBLY: "$set
echo "--------------------------------------------------------------------------------"
mkdir -p "$assembly_dir"
# Run Porechop
#   --no_split: to save time, chimeras shouldn't matter for Rebaler anyway
#   --check_reads 1000: to save time, it's all one barcode, so 1000 reads should be plenty
porechop -i $1 -o $trimmed_reads --no_split --threads $threads --check_reads 1000
if [ "$read_type" = "fasta" ]; then
    seqtk seq $trimmed_reads > "$assembly_dir"/temp.fasta  # make one-line-per-seq fasta
    mv "$assembly_dir"/temp.fasta $trimmed_reads
fi
# Run multiple replicates of Racon on shuffled reads and a rotated reference.
# Each assembly is then coursely shredded into 'reads' for the final Racon round.
for i in {0..9}; do
    sample_reads="$assembly_dir"/"$set"_sample_"$i"."$read_type"
    sample_assembly="$assembly_dir"/"$set"_assembly_"$i".fasta
    sample_reference="$assembly_dir"/"$set"_reference_"$i".fasta
    if [ "$read_type" = "fastq" ]; then
        cat $trimmed_reads | paste - - - - | shuf | awk -F'\t' '{OFS="\n"; print $1,$2,$3,$4}' > $sample_reads
    elif [ "$read_type" = "fasta" ]; then
        cat $trimmed_reads | paste - - | shuf | awk -F'\t' '{OFS="\n"; print $1,$2}' > $sample_reads
    fi
    "$scripts_dir"/rotate_reference.py reference.fasta $i > $sample_reference
    rebaler -t $threads $sample_reference $sample_reads > $sample_assembly
    "$scripts_dir"/shred_assembly.py $sample_assembly $i 50000 >> $assembly_reads
    rm $sample_reads $sample_assembly $sample_reference
done
rm $trimmed_reads
rebaler -t $threads reference.fasta $assembly_reads > $final_assembly
rm $assembly_reads




printf "\n\n\n\n"
echo "ASSESS ASSEMBLY: "$set
echo "--------------------------------------------------------------------------------"
python3 "$scripts_dir"/chop_up_assembly.py $final_assembly 10000 > $assembly_pieces
minimap2 -x asm5 -t $threads -c reference.fasta $assembly_pieces > $assembly_alignment
pypy3 "$scripts_dir"/read_length_identity.py $assembly_pieces $assembly_alignment > $assembly_data
rm $assembly_pieces $assembly_alignment
if [ ! -f "$results_dir"/assembly_error_details ]; then
    printf "assembly\tdcm\thomo_del\thomo_ins\tother_del\tother_ins\tsub\n" > "$results_dir"/assembly_error_details
fi
printf $set"\t" >> "$results_dir"/assembly_error_details
nucmer --prefix="$prefix" reference.fasta $final_assembly
delta-filter -r -q "$prefix".delta > "$prefix".filter
# show-snps -ClrTH -x5 "$prefix".filter | python3 "$scripts_dir"/error_summary.py "$ref_contig" "$assembly_contig" >> "$results_dir"/assembly_error_details
printf $set"\tassembly\t" >> "$results_dir"/substitution_counts
show-snps -ClrTH "$prefix".filter | awk '$2 != "." && $3 != "."' | wc -l >> "$results_dir"/substitution_counts
rm "$prefix".delta "$prefix".filter




printf "\n\n\n\n"
echo "NANOPOLISH: "$set
echo "--------------------------------------------------------------------------------"
mkdir -p $nanopolish_dir
if [ -f $sequencing_summary ]; then
    echo "Running nanopolish index with sequencing_summary.txt file"
    nanopolish index -d fast5s -s $sequencing_summary $1
else
    echo "Running nanopolish index without sequencing_summary.txt file"
    nanopolish index -d fast5s $1
fi
minimap2 -x map-ont -a -t $threads $final_assembly $1 | samtools sort > $bam_file
samtools index $bam_file
nanopolish_makerange.py $final_assembly | parallel --results "$nanopolish_dir"/"$set"_nanopolish.results -P $threads nanopolish variants --consensus -o "$nanopolish_dir"/"$set"_polished.{1}.vcf -w {1} -r $1 -b $bam_file -g $final_assembly -t 2 --min-candidate-frequency 0.1 --methylation-aware=dcm,dam --fix-homopolymers
nanopolish vcf2fasta -g $final_assembly "$nanopolish_dir"/"$set"_polished.*.vcf > $nanopolish_assembly
rm "$nanopolish_dir"/"$set"_polished.*.vcf "$nanopolish_dir"/"$set"_nanopolish.results*
rm "$1".index "$1".index.fai "$1".index.gzi "$1".index.readdb $bam_file "$bam_file".bai
rm "$final_assembly".fai




printf "\n\n\n\n"
echo "ASSESS NANOPOLISHED ASSEMBLY: "$set
echo "--------------------------------------------------------------------------------"
python3 "$scripts_dir"/chop_up_assembly.py $nanopolish_assembly 10000 > $nanopolish_pieces
minimap2 -x asm5 -t $threads -c reference.fasta $nanopolish_pieces > $nanopolish_alignment
pypy3 "$scripts_dir"/read_length_identity.py $nanopolish_pieces $nanopolish_alignment > $nanopolish_data
rm $nanopolish_pieces $nanopolish_alignment
if [ ! -f "$results_dir"/nanopolish_error_details ]; then
    printf "assembly\tdcm\thomo_del\thomo_ins\tother_del\tother_ins\tsub\n" > "$results_dir"/nanopolish_error_details
fi
printf $set"\t" >> "$results_dir"/nanopolish_error_details
nucmer --prefix="$prefix" reference.fasta $nanopolish_assembly
delta-filter -r -q "$prefix".delta > "$prefix".filter
show-snps -ClrTH -x5 "$prefix".filter | python3 "$scripts_dir"/error_summary.py "$ref_contig" "$assembly_contig" >> "$results_dir"/nanopolish_error_details
printf $set"\tnanopolish\t" >> "$results_dir"/substitution_counts
show-snps -ClrTH "$prefix".filter | awk '$2 != "." && $3 != "."' | wc -l >> "$results_dir"/substitution_counts
rm "$prefix".delta "$prefix".filter




printf "\n\n\n\n"
echo "NANOPOLISH (SECOND TIME): "$set
echo "--------------------------------------------------------------------------------"
if [ -f $sequencing_summary ]; then
    echo "Running nanopolish index with sequencing_summary.txt file"
    nanopolish index -d fast5s -s $sequencing_summary $1
else
    echo "Running nanopolish index without sequencing_summary.txt file"
    nanopolish index -d fast5s $1
fi
minimap2 -x map-ont -a -t $threads $nanopolish_assembly $1 | samtools sort > $bam_file_2
samtools index $bam_file_2
nanopolish_makerange.py $nanopolish_assembly | parallel --results "$nanopolish_dir"/"$set"_2_nanopolish.results -P $threads nanopolish variants --consensus -o "$nanopolish_dir"/"$set"_2_polished.{1}.vcf -w {1} -r $1 -b $bam_file_2 -g $nanopolish_assembly -t 2 --min-candidate-frequency 0.1 --methylation-aware=dcm,dam --fix-homopolymers
nanopolish vcf2fasta -g $nanopolish_assembly "$nanopolish_dir"/"$set"_2_polished.*.vcf > $nanopolish_2_assembly
rm "$nanopolish_dir"/"$set"_2_polished.*.vcf "$nanopolish_dir"/"$set"_2_nanopolish.results*
rm "$1".index "$1".index.fai "$1".index.gzi "$1".index.readdb $bam_file_2 "$bam_file_2".bai
rm "$nanopolish_assembly".fai



printf "\n\n\n\n"
echo "ASSESS NANOPOLISHED (SECOND TIME) ASSEMBLY: "$set
echo "--------------------------------------------------------------------------------"
python3 "$scripts_dir"/chop_up_assembly.py $nanopolish_2_assembly 10000 > $nanopolish_2_pieces
minimap2 -x asm5 -t $threads -c reference.fasta $nanopolish_2_pieces > $nanopolish_2_alignment
pypy3 "$scripts_dir"/read_length_identity.py $nanopolish_2_pieces $nanopolish_2_alignment > $nanopolish_2_data
rm $nanopolish_2_pieces $nanopolish_2_alignment
if [ ! -f "$results_dir"/nanopolish_2_error_details ]; then
    printf "assembly\tdcm\thomo_del\thomo_ins\tother_del\tother_ins\tsub\n" > "$results_dir"/nanopolish_2_error_details
fi
printf $set"\t" >> "$results_dir"/nanopolish_2_error_details
nucmer --prefix="$prefix"_2 reference.fasta $nanopolish_2_assembly
delta-filter -r -q "$prefix"_2.delta > "$prefix"_2.filter
show-snps -ClrTH -x5 "$prefix"_2.filter | python3 "$scripts_dir"/error_summary.py "$ref_contig" "$assembly_contig" >> "$results_dir"/nanopolish_2_error_details
rm "$prefix"_2.delta "$prefix"_2.filter




printf "\n\n\n\n"
echo "NANOPOLISH (THIRD TIME): "$set
echo "--------------------------------------------------------------------------------"
if [ -f $sequencing_summary ]; then
    echo "Running nanopolish index with sequencing_summary.txt file"
    nanopolish index -d fast5s -s $sequencing_summary $1
else
    echo "Running nanopolish index without sequencing_summary.txt file"
    nanopolish index -d fast5s $1
fi
minimap2 -x map-ont -a -t $threads $nanopolish_2_assembly $1 | samtools sort > $bam_file_3
samtools index $bam_file_3
nanopolish_makerange.py $nanopolish_2_assembly | parallel --results "$nanopolish_dir"/"$set"_3_nanopolish.results -P $threads nanopolish variants --consensus -o "$nanopolish_dir"/"$set"_3_polished.{1}.vcf -w {1} -r $1 -b $bam_file_3 -g $nanopolish_2_assembly -t 2 --min-candidate-frequency 0.1 --methylation-aware=dcm,dam --fix-homopolymers
nanopolish vcf2fasta -g $nanopolish_2_assembly "$nanopolish_dir"/"$set"_3_polished.*.vcf > $nanopolish_3_assembly
rm "$nanopolish_dir"/"$set"_3_polished.*.vcf "$nanopolish_dir"/"$set"_3_nanopolish.results*
rm "$1".index "$1".index.fai "$1".index.gzi "$1".index.readdb $bam_file_3 "$bam_file_3".bai
rm "$nanopolish_2_assembly".fai



printf "\n\n\n\n"
echo "ASSESS NANOPOLISHED (THIRD TIME) ASSEMBLY: "$set
echo "--------------------------------------------------------------------------------"
python3 "$scripts_dir"/chop_up_assembly.py $nanopolish_3_assembly 10000 > $nanopolish_3_pieces
minimap2 -x asm5 -t $threads -c reference.fasta $nanopolish_3_pieces > $nanopolish_3_alignment
pypy3 "$scripts_dir"/read_length_identity.py $nanopolish_3_pieces $nanopolish_3_alignment > $nanopolish_3_data
rm $nanopolish_3_pieces $nanopolish_3_alignment
if [ ! -f "$results_dir"/nanopolish_3_error_details ]; then
    printf "assembly\tdcm\thomo_del\thomo_ins\tother_del\tother_ins\tsub\n" > "$results_dir"/nanopolish_3_error_details
fi
printf $set"\t" >> "$results_dir"/nanopolish_3_error_details
nucmer --prefix="$prefix"_3 reference.fasta $nanopolish_3_assembly
delta-filter -r -q "$prefix"_3.delta > "$prefix"_3.filter
show-snps -ClrTH -x5 "$prefix"_3.filter | python3 "$scripts_dir"/error_summary.py "$ref_contig" "$assembly_contig" >> "$results_dir"/nanopolish_3_error_details
rm "$prefix"_3.delta "$prefix"_3.filter




printf "\n\n\n\n"
echo "NANOPOLISH (FOURTH TIME): "$set
echo "--------------------------------------------------------------------------------"
if [ -f $sequencing_summary ]; then
    echo "Running nanopolish index with sequencing_summary.txt file"
    nanopolish index -d fast5s -s $sequencing_summary $1
else
    echo "Running nanopolish index without sequencing_summary.txt file"
    nanopolish index -d fast5s $1
fi
minimap2 -x map-ont -a -t $threads $nanopolish_3_assembly $1 | samtools sort > $bam_file_4
samtools index $bam_file_4
nanopolish_makerange.py $nanopolish_3_assembly | parallel --results "$nanopolish_dir"/"$set"_4_nanopolish.results -P $threads nanopolish variants --consensus -o "$nanopolish_dir"/"$set"_4_polished.{1}.vcf -w {1} -r $1 -b $bam_file_4 -g $nanopolish_3_assembly -t 2 --min-candidate-frequency 0.1 --methylation-aware=dcm,dam --fix-homopolymers
nanopolish vcf2fasta -g $nanopolish_3_assembly "$nanopolish_dir"/"$set"_4_polished.*.vcf > $nanopolish_4_assembly
rm "$nanopolish_dir"/"$set"_4_polished.*.vcf "$nanopolish_dir"/"$set"_4_nanopolish.results*
rm "$1".index "$1".index.fai "$1".index.gzi "$1".index.readdb $bam_file_4 "$bam_file_4".bai
rm "$nanopolish_3_assembly".fai



printf "\n\n\n\n"
echo "ASSESS NANOPOLISHED (FOURTH TIME) ASSEMBLY: "$set
echo "--------------------------------------------------------------------------------"
python3 "$scripts_dir"/chop_up_assembly.py $nanopolish_4_assembly 10000 > $nanopolish_4_pieces
minimap2 -x asm5 -t $threads -c reference.fasta $nanopolish_4_pieces > $nanopolish_4_alignment
pypy3 "$scripts_dir"/read_length_identity.py $nanopolish_4_pieces $nanopolish_4_alignment > $nanopolish_4_data
rm $nanopolish_4_pieces $nanopolish_4_alignment
if [ ! -f "$results_dir"/nanopolish_4_error_details ]; then
    printf "assembly\tdcm\thomo_del\thomo_ins\tother_del\tother_ins\tsub\n" > "$results_dir"/nanopolish_4_error_details
fi
printf $set"\t" >> "$results_dir"/nanopolish_4_error_details
nucmer --prefix="$prefix"_4 reference.fasta $nanopolish_4_assembly
delta-filter -r -q "$prefix"_4.delta > "$prefix"_4.filter
show-snps -ClrTH -x5 "$prefix"_4.filter | python3 "$scripts_dir"/error_summary.py "$ref_contig" "$assembly_contig" >> "$results_dir"/nanopolish_4_error_details
rm "$prefix"_4.delta "$prefix"_4.filter
