fast5_dir=/var/lib/MinKNOW/data/basecalling_comparison/fast5s
output_prefix=/var/lib/MinKNOW/data/basecalling_comparison/basecalling_output
guppy_installers=/var/lib/MinKNOW/data/basecalling_comparison/guppy_installers

cd $output_prefix

# Just using the latest version of Guppy for the custom models.
v=2.2.3

apt install "$guppy_installers"/ont-guppy_"$v"-1~xenial_amd64.deb

out_dir="$output_prefix"/guppy_v"$v"_custom_kp_model
/usr/bin/time -v -o guppy_v"$v"_custom_kp_model.time guppy_basecaller --model_file holtlab_kp_r9.4_r9.4.1_nov_2018.jsn --config dna_r9.4.1_450bps.cfg --device auto -i $fast5_dir -t 12 -s $out_dir
cat "$out_dir"/*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > guppy_v"$v"_custom_kp_model.fastq

out_dir="$output_prefix"/guppy_v"$v"_custom_kp_big_net_model
/usr/bin/time -v -o guppy_v"$v"_custom_kp_big_net_model.time guppy_basecaller --model_file holtlab_kp_big_r9.4_r9.4.1_nov_2018.jsn --config dna_r9.4.1_450bps.cfg -i $fast5_dir -t 12 -s $out_dir
cat "$out_dir"/*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > guppy_v"$v"_custom_kp_big_net_model.fastq
