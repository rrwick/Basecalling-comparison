fast5_dir=/var/lib/MinKNOW/data/basecalling_comparison/fast5s
output_prefix=/var/lib/MinKNOW/data/basecalling_comparison/basecalling_output
albacore_installers=/var/lib/MinKNOW/data/basecalling_comparison/albacore_installers

cd $output_prefix


for v in 2.1.0 2.1.1 2.1.1 2.1.2 2.1.3 2.1.5 2.1.6 2.1.7 2.1.9 2.1.10 2.2.1 2.2.2 2.2.4 2.2.6 2.2.7 2.3.0 2.3.1 2.3.2 2.3.3 2.3.4; do
    pip3 install "$albacore_installers"/ont_albacore-"$v"-cp35-cp35m-manylinux1_x86_64.whl
    out_dir="$output_prefix"/albacore_v"$v"
    /usr/bin/time -v -o albacore_v"$v".time read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -i $fast5_dir -t 12 -s $out_dir -o fastq --disable_filtering --disable_pings
    cat "$out_dir"/workspace/*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > albacore_v"$v".fastq
done

# These versions don't have --disable_pings.
for v in 2.0.0 2.0.1 2.0.2; do
    pip3 install "$albacore_installers"/ont_albacore-"$v"-cp35-cp35m-manylinux1_x86_64.whl
    out_dir="$output_prefix"/albacore_v"$v"
    /usr/bin/time -v -o albacore_v"$v".time read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -i $fast5_dir -t 12 -s $out_dir -o fastq --disable_filtering
    cat "$out_dir"/workspace/*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > albacore_v"$v".fastq
done

# These versions don't have --disable_filtering.
for v in 1.2.0 1.2.1 1.2.2 1.2.3 1.2.4 1.2.5 1.2.6 1.3.1; do
    pip3 install "$albacore_installers"/ont_albacore-"$v"-cp35-cp35m-manylinux1_x86_64.whl
    out_dir="$output_prefix"/albacore_v"$v"
    /usr/bin/time -v -o albacore_v"$v".time read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -i $fast5_dir -t 12 -s $out_dir -o fastq
    cat "$out_dir"/workspace/*.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > albacore_v"$v".fastq
done

# These versions mess up the fastq files (extra line break).
for v in 1.1.0 1.1.1 1.1.2; do
    pip3 install "$albacore_installers"/ont_albacore-"$v"-cp35-cp35m-manylinux1_x86_64.whl
    out_dir="$output_prefix"/albacore_v"$v"
    /usr/bin/time -v -o albacore_v"$v".time read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -i $fast5_dir -t 12 -s $out_dir -o fastq
    cat "$out_dir"/workspace/*.fastq | paste - - - - - | sort -k1,1 -t " " | tr "\t" "\n" > albacore_v"$v".fastq
done

# These versions can't output to fastq.
for v in 1.0.1 1.0.2 1.0.3 1.0.4; do
    pip3 install "$albacore_installers"/ont_albacore-"$v"-cp35-cp35m-manylinux1_x86_64.whl
    out_dir="$output_prefix"/albacore_v"$v"
    /usr/bin/time -v -o albacore_v"$v".time read_fast5_basecaller.py -f FLO-MIN106 -k SQK-LSK108 -i $fast5_dir -t 12 -s $out_dir
    ~/Fast5-to-Fastq/fast5_to_fastq.py "$out_dir"/workspace | paste - - - - | sed 's/_Basecall_1D_template/ Basecall_1D_template/' | sort -k1,1 -t " " | tr "\t" "\n" > albacore_v"$v".fastq
done

# These versions don't have flowcell and kit options (need to specify config).
for v in 0.8.2 0.8.4 0.9.1; do
    pip3 install "$albacore_installers"/ont_albacore-"$v"-cp35-cp35m-manylinux1_x86_64.whl
    out_dir="$output_prefix"/albacore_v"$v"
    /usr/bin/time -v -o albacore_v"$v".time read_fast5_basecaller.py -c FLO-MIN106_LSK108_linear.cfg -i $fast5_dir -t 12 -s $out_dir
    ~/Fast5-to-Fastq/fast5_to_fastq.py "$out_dir"/workspace | paste - - - - | sed 's/_Basecall_1D_template/ Basecall_1D_template/' | sort -k1,1 -t " " | tr "\t" "\n" > albacore_v"$v".fastq
done

# v0.7.5 is the only one without a Python wheel.
v=0.7.5
apt install "$albacore_installers"/python3-ont-albacore_"$v"-1-xenial_amd64.deb
out_dir="$output_prefix"/albacore_v"$v"
/usr/bin/time -v -o albacore_v"$v".time read_fast5_basecaller.py -c FLO-MIN106_LSK108_linear.cfg -i $fast5_dir -t 12 -s $out_dir
~/Fast5-to-Fastq/fast5_to_fastq.py "$out_dir"/workspace | paste - - - - | sed 's/_Basecall_1D_template/ Basecall_1D_template/' | sort -k1,1 -t " " | tr "\t" "\n" > albacore_v"$v".fastq
