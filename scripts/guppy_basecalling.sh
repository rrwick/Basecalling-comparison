fast5_dir=/var/lib/MinKNOW/data/basecalling_comparison/fast5s
output_prefix=/var/lib/MinKNOW/data/basecalling_comparison/basecalling_output
guppy_installers=/var/lib/MinKNOW/data/basecalling_comparison/guppy_installers

cd $output_prefix


# For the most recent versions of Guppy, I'm using both the default and flip-flop models.
for v in 2.1.3 2.2.3; do
    apt install "$guppy_installers"/ont-guppy_"$v"-1~xenial_amd64.deb
    out_dir="$output_prefix"/guppy_v"$v"
    /usr/bin/time -v -o guppy_v"$v".time guppy_basecaller --config dna_r9.4.1_450bps.cfg --device auto -i $fast5_dir -t 12 -s $out_dir
    cat "$out_dir"/*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > guppy_v"$v".fastq
    out_dir="$output_prefix"/guppy_v"$v"_flipflop
    /usr/bin/time -v -o guppy_v"$v"_flipflop.time guppy_basecaller --config dna_r9.4.1_450bps_flipflop.cfg --device auto -i $fast5_dir -t 12 -s $out_dir
    cat "$out_dir"/*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > guppy_v"$v"_flipflop.fastq
done

for v in 1.5.1 1.6.0 1.8.1 1.8.3 1.8.5; do
    apt install "$guppy_installers"/ont-guppy_"$v"-1~xenial_amd64.deb
    out_dir="$output_prefix"/guppy_v"$v"
    /usr/bin/time -v -o guppy_v"$v".time guppy_basecaller --config dna_r9.4_450bps.cfg --device auto -i $fast5_dir -t 12 -s $out_dir
    cat "$out_dir"/*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > guppy_v"$v".fastq
done

# I couldn't make v1.4.3 run on the GPU so I used the CPU instead.
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

# I couldn't make v0.3.0 run on the GPU so I used the CPU instead.
v=0.3.0
apt install "$guppy_installers"/ont-guppy_"$v"-1~xenial_amd64.deb
out_dir="$output_prefix"/guppy_v"$v"
/usr/bin/time -v -o guppy_v"$v".time guppy --config dna_r9.4_450bps.cfg --device cpu -i $fast5_dir -s $out_dir
cat "$out_dir"/*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > guppy_v"$v".fastq
