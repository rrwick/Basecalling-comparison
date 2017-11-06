#!/usr/bin/env python3
"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Basecalling-comparison

This script is a Nanopolish wrapper I wrote for use on my SLURM-managed cluster. It does the read
alignment, launches Nanopolish jobs, waits for them to finish and merges them together. If any
parts of the assembly fail in Nanopolish it replaces them with Ns so the merge can complete.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should
have received a copy of the GNU General Public License along with this program. If not, see
<http://www.gnu.org/licenses/>.
"""

import sys
import os
import subprocess
import time
import shutil


def main():
    assembly_filename = os.path.abspath(sys.argv[1])
    read_filename = os.path.abspath(sys.argv[2])
    raw_fast5_dir = os.path.abspath(sys.argv[3])
    output_dir = os.path.abspath(sys.argv[4])
    nanopolish_dir = os.path.abspath(sys.argv[5])
    threads = int(sys.argv[6])

    methylation_aware = False
    try:
        if sys.argv[7].lower() == 'meth':
            methylation_aware = True
    except IndexError:
        pass

    nanopolish_exec = os.path.join(nanopolish_dir, 'nanopolish')
    nanopolish_makerange = os.path.join(nanopolish_dir, 'scripts', 'nanopolish_makerange.py')
    nanopolish_merge = os.path.join(nanopolish_dir, 'scripts', 'nanopolish_merge.py')

    set_name = assembly_filename.split('/')[-1].split('.fasta')[0]
    print('\nPreparing to run Nanopolish for ' + set_name)

    final_assembly = ('../' + set_name + '.fasta').replace('_assembly.fasta', '_nanopolish.fasta')
    if methylation_aware:
        final_assembly = final_assembly.replace('_nanopolish.fasta', '_nanopolish_meth.fasta')

    pid = str(os.getpid())
    temp_dir = os.path.join(output_dir, pid + '_temp_dir')
    os.mkdir(temp_dir)
    print('Moving into ' + temp_dir)
    os.chdir(temp_dir)

    print('Getting ranges: ', end='')
    polish_ranges = get_nanopolish_ranges(nanopolish_makerange, assembly_filename)
    print(', '.join(polish_ranges))

    # # Align reads with BWA MEM
    # print('Aligning reads')
    # index_command = 'bwa index ' + assembly_filename
    # subprocess.run(index_command, shell=True, check=True)
    # alignment_command = 'bwa mem -x ont2d -t ' + str(threads) + ' ' + assembly_filename + ' ' + read_filename + ' | samtools sort -o reads.sorted.bam -T reads.tmp -'
    # subprocess.run(alignment_command, shell=True, check=True)
    # subprocess.run('samtools index reads.sorted.bam', shell=True, check=True)

    # Align reads with minimap2
    print('Aligning reads')
    alignment_command = ('minimap2 -x map10k -a -t ' + str(threads) + ' ' + 
                         assembly_filename + ' ' + read_filename + 
                         ' | samtools sort -o reads.sorted.bam -T reads.tmp -')
    subprocess.run(alignment_command, shell=True, check=True)
    subprocess.run('samtools index reads.sorted.bam', shell=True, check=True)

    # Run Nanopolish index on reads
    print('Running nanopolish index:')
    index_command = nanopolish_exec + ' index -d ' + raw_fast5_dir + ' ' + read_filename
    subprocess.run(index_command, shell=True, check=True)

    # Run Nanopolish variants on ranges
    print('Launching SLURM jobs:', flush=True)
    job_prefix = 'Nanopolish_' + pid + '_'
    for polish_range in polish_ranges:
        job_name = job_prefix + polish_range
        variants_command = nanopolish_exec + ' variants --consensus polished.' + polish_range + '.fa -w ' + polish_range + ' -r ' + read_filename + ' -b reads.sorted.bam -g ' + assembly_filename + ' -t 2 --min-candidate-frequency 0.1'
        if methylation_aware:
            variants_command += ' --methylation-aware=dcm,dam'
        sbatch_command = 'sbatch -p sysgen --nodes=1 --job-name=' + job_name + ' --ntasks=1 --cpus-per-task=2 --mem=4096 --time=0-4:0:00 --wrap "' + variants_command + '"'
        print(sbatch_command)
        subprocess.run(sbatch_command, shell=True, check=True)

    # Wait for jobs to finish
    start_time = time.time()
    while True:
        time.sleep(60)
        remaining_jobs = get_remaining_nanopolish_jobs(job_prefix)
        if remaining_jobs == 0:
            print('All Nanopolish jobs are done!')
            break
        elapsed_time = str(int(round(time.time() - start_time)))
        print('Waiting for Nanopolish jobs to finish... (' + elapsed_time + ' sec elapsed, ' + str(remaining_jobs) + ' jobs remaining)', flush=True)

    # Make empty ranges, if necessary
    incomplete_ranges = [x for x in polish_ranges if not os.path.isfile('polished.' + x + '.fa')]
    if incomplete_ranges:
        print('WARNING: some ranges did not complete: ' + ', '.join(incomplete_ranges))
    for incomplete_range in incomplete_ranges:
        fasta_filename = 'polished.' + incomplete_range + '.fa'
        start = int(incomplete_range.split(':')[-1].split('-')[0])
        end = int(incomplete_range.split('-')[-1])
        range_size = end - start
        with open(fasta_filename, 'wt') as fasta:
            fasta.write('>')
            fasta.write(incomplete_range)
            fasta.write('\n')
            fasta.write('N' * range_size)
            fasta.write('\n')

    # Merge results together
    merge_command = 'python ' + nanopolish_merge + ' polished.*.fa > ' + final_assembly
    subprocess.run(merge_command, shell=True, check=True)

    os.chdir('..')
    shutil.rmtree(temp_dir)


def get_nanopolish_ranges(nanopolish_makerange, assembly_filename):
    command = 'python ' + nanopolish_makerange + ' ' + assembly_filename
    range_out = subprocess.check_output(command, shell=True).splitlines()
    return [x.decode() for x in range_out]


def get_remaining_nanopolish_jobs(job_prefix):
    current_jobs = subprocess.check_output('squeue -o "%.70j %.8i %.10T"', shell=True).decode().splitlines()
    remaining_jobs = [x for x in current_jobs if job_prefix in x]
    return len(remaining_jobs)


if __name__ == '__main__':
    main()
