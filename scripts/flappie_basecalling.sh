fast5_dir=/var/lib/MinKNOW/data/basecalling_comparison/fast5s
output_prefix=/var/lib/MinKNOW/data/basecalling_comparison/basecalling_output

# I built all of the Scrappie versions and copied their binaries to a single
# directory. Most were easy to build, but versions prior to v0.2.6 had an
# issue with finding the hdf5.h header. I had to copy Tim's CMakeLists.txt fix
# (in v0.2.6) into those versions to make them build. Also versions 0.0.2 and
# 0.0.3 didn't use CMake, so I needed to make a few tweaks to their Makefiles.
flappie_dir=/home/holt-ont/Flappie

cd $output_prefix

export OPENBLAS_NUM_THREADS=1

for v in 1.0 1.1; do
    out_dir="$output_prefix"/flappie_v"$v"
    flappie="$flappie_dir"/flappie-"$v"/flappie
    # I can't use time for Flappie because I run it in separate processes with the parallel
    # command, so I just run date before and after to get the total time.
    date
    find $fast5_dir -name \*.fast5 | parallel -P 12 -X $flappie --model r941_native > flappie_v"$v".fastq
    date
done