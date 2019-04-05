#!/usr/bin/env python3
"""
Copyright 2019 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Basecalling-comparison

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should
have received a copy of the GNU General Public License along with this program. If not, see
<http://www.gnu.org/licenses/>.

This script moves fast5 files into subdirectories with a particular number of files per directory.
"""

import glob
import math
import os
import shutil
import sys

input_dir = sys.argv[1]
target_reads_per_dir = int(sys.argv[2])

fast5_files = glob.glob(input_dir + '/*.fast5')
file_count = len(fast5_files)
print('Found {} fast5 files in {}'.format(file_count, input_dir))

dir_count = int(round(file_count / target_reads_per_dir))
if dir_count == 0:
    dir_count = 1
files_per_dir = int(math.ceil(file_count / dir_count))

dir_num = 0
while file_count > 0:
    dir_name = input_dir + '/{:04d}'.format(dir_num)
    os.makedirs(dir_name)

    if file_count >= files_per_dir:
        move_files = files_per_dir
    else:
        move_files = file_count
    print('Moving {} files to {}'.format(move_files, dir_name))
    for f in fast5_files[:move_files]:
        shutil.move(f, dir_name)

    fast5_files = fast5_files[move_files:]
    file_count -= move_files
    dir_num += 1
