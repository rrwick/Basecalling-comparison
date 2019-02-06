fast5_dir=/var/lib/MinKNOW/data/basecalling_comparison/fast5s
output_prefix=/var/lib/MinKNOW/data/basecalling_comparison/basecalling_output

# I built all of the Scrappie versions and copied their binaries to a single directory:
binaries=/home/holt-ont/Scrappie/binaries

cd $output_prefix

export OMP_NUM_THREADS=$(nproc)
export OPENBLAS_NUM_THREADS=1


for v in 1.4.1 1.4.0 1.3.2 1.3.1 1.3.0 1.2.0 1.1.1 1.1.0; do
    out_dir="$output_prefix"/scrappie_v"$v"
    scrappie="$binaries"/scrappie-v"$v"
    /usr/bin/time -v -o scrappie_raw_v"$v".time $scrappie raw $fast5_dir > scrappie_raw_v"$v".fasta
    /usr/bin/time -v -o scrappie_events_v"$v".time $scrappie events $fast5_dir > scrappie_events_v"$v".fasta
done


# These versions can't call events on their own - they rely on events already
# being present, e.g. from Albacore. I therefore skipped scrappie events
# (because it's not a stand-alone basecaller) and just ran scrappie raw.
for v in 1.0.0 0.3.2 0.3.1 0.3.0; do
    out_dir="$output_prefix"/scrappie_v"$v"
    scrappie="$binaries"/scrappie-v"$v"
    /usr/bin/time -v -o scrappie_raw_v"$v".time $scrappie raw $fast5_dir > scrappie_raw_v"$v".fasta
done


# Versions earlier than 0.3.0 can only do event-based basecalling, so I'm skipping all of them.
