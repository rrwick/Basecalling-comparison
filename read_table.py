#!/usr/bin/env python3
"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Basecalling-comparison

This script takes a single argument (a directory containing fast5 files) and prints a table to
stdout with some basic read info.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should
have received a copy of the GNU General Public License along with this program. If not, see
<http://www.gnu.org/licenses/>.
"""

import h5py
import re
import sys
import os


def main():
    print('\t'.join(['Name', 'Fast5 name', 'Run name', 'Signal length', 'Start time']))
    rows = []
    count = 0
    for fast5_filename in find_all_fast5s(sys.argv[1]):
        fast5_name = fast5_filename.split('/')[-1][:-6]
        with h5py.File(fast5_filename) as read:
            names = get_hdf5_names(read)
            read_id = get_read_id(read, names)
            run_name = get_run_name(read, names)
            signal_length = get_signal_length(read, names)
            start_time = get_start_time(read, names)
            count += 1
            if count % 10 == 0:
                print('.', file=sys.stderr, flush=True, end='')
            rows.append([read_id, fast5_name, run_name, str(signal_length), str(start_time)])
    rows = sorted(rows)
    print('', file=sys.stderr, flush=True)
    for row in rows:
        print('\t'.join(row))


def find_all_fast5s(directory):
    fast5s = []
    for dir_name, _, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith('.fast5'):
                fast5s.append(os.path.join(dir_name, filename))
    return fast5s


def get_read_id(fast5_file, names):
    r = re.compile(r'/Read_\d+$')
    info_locations = [x for x in names if r.search(x)]
    assert len(info_locations) == 1
    info_location = info_locations[0]
    assert 'read_id' in fast5_file[info_location].attrs
    return fast5_file[info_location].attrs['read_id'].decode()


def get_run_name(fast5_file, names):
    tracking_names = [x for x in names if x.endswith('tracking_id')]
    assert len(tracking_names) == 1
    tracking_location = tracking_names[0]
    assert 'sample_id' in fast5_file[tracking_location].attrs
    assert 'exp_script_purpose' in fast5_file[tracking_location].attrs
    sample_id = fast5_file[tracking_location].attrs['sample_id'].decode()
    exp_script_purpose = fast5_file[tracking_location].attrs['exp_script_purpose'].decode()
    return sample_id + '_' + exp_script_purpose


def get_signal_length(fast5_file, names):
    signal_names = [x for x in names if x.endswith('Signal')]
    assert len(signal_names) == 1
    signal_location = signal_names[0]
    return fast5_file[signal_location].shape[0]


def get_start_time(fast5_file, names):
    r = re.compile(r'/Read_\d+$')
    info_locations = [x for x in names if r.search(x)]
    assert len(info_locations) == 1
    info_location = info_locations[0]
    assert 'start_time' in fast5_file[info_location].attrs
    return fast5_file[info_location].attrs['start_time']


def get_hdf5_names(fast5_file):
    names = []
    fast5_file.visit(names.append)
    return names


if __name__ == '__main__':
    main()
