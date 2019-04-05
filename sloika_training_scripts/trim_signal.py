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

This is a tool for trimming fast5s at the signal level.

It does two steps for each read:
  1) Trim off what looks like open-pore signal from the start and end of reads.
  2) Trim off a fixed amount of signal from the start and end of reads.

The final result should be reads with no adapter sequence, making them more suitable for Sloika.
"""

import argparse
import glob
import h5py
import numpy as np
import os
import shutil
import sys


def get_arguments():
    parser = argparse.ArgumentParser(description='Trim fast5 files at the signal level')

    parser.add_argument('--trim_amount', type=int, required=False, default=2000,
                        help='Signal values to be trimmed from the start/end of the signal')
    parser.add_argument('--min_size', type=int, required=False, default=50000,
                        help='Don\'t keep reads with a trimmed signal smaller than this')

    parser.add_argument('in_dir', type=str,
                        help='Directory containing fast5 files to trim')
    parser.add_argument('out_dir', type=str,
                        help='Output directory of trimmed fast5 files (will be made)')

    args = parser.parse_args()
    return args


def main():
    args = get_arguments()
    make_output_dir(args.out_dir)

    for original_file in glob.iglob(args.in_dir + '/**/*.fast5', recursive=True):
        print(original_file)

        # Make a copy of the file which we will modify.       
        new_file = '{}/{}'.format(args.out_dir, original_file.rpartition('/')[2])
        print('    making a copy here: {}'.format(new_file))
        if os.path.exists(new_file):
            print()
            sys.exit('Error: {} already exists'.format(new_file))
        shutil.copyfile(original_file, new_file)

        # Find the signal.
        hdf5_file = h5py.File(new_file, 'r+')
        signal_location = get_signal_location(hdf5_file)
        read_location = signal_location.replace('/Signal', '')
        signal = hdf5_file[signal_location].value
        orignal_length = len(signal)
        print('    original signal length = {:,}'.format(orignal_length))

        # Trim off open-pore signal
        try:
            start_trim = find_signal_start_pos(signal)
            print('    open-pore signal at start: {}'.format(start_trim))
            end_trim = find_signal_start_pos(signal[::-1])
            print('    open-pore signal at end: {}'.format(end_trim))
        except CannotTrim:
            print('    cannot trim - skipping')
            hdf5_file.close()
            os.remove(new_file)
            continue

        if start_trim + end_trim >= orignal_length:
            print('    too short - skipping')
            hdf5_file.close()
            os.remove(new_file)
            continue
        signal = signal[start_trim:-end_trim]
        no_open_pore_length = len(signal)
        print('    open-pore trimmed signal length = {:,}'.format(no_open_pore_length))

        if no_open_pore_length - (2 * args.trim_amount) < args.min_size:
            print('    too short - skipping')
            hdf5_file.close()
            os.remove(new_file)
        else:
            # Trim the signal and rewrite it.
            signal = signal[args.trim_amount:-args.trim_amount]
            del hdf5_file[signal_location]
            dset = hdf5_file.create_dataset(signal_location, compression='gzip',
                                            compression_opts=9, data=signal)
            new_length = len(signal)
            hdf5_file[read_location].attrs['duration'] = new_length
            print('    final trimmed signal length = {:,}'.format(new_length))
            hdf5_file.close()
        print()


def make_output_dir(out_dir):
    if os.path.exists(out_dir):
        sys.exit('Error: output directory already exists')
    os.makedirs(out_dir)


def get_signal_location(hdf5_file):
    names = []
    hdf5_file.visit(names.append)
    signal_locations = sorted([x for x in names if x.endswith('/Signal')])
    assert len(signal_locations) == 1
    return signal_locations[0]


def find_signal_start_pos(signal):
    """
    Given a signal, this function attempts to identify the approximate position where the open
    pore signal ends and the real signal begins.
    """
    initial_trim_size = 10
    trim_increment = 25
    stdev_threshold = 20
    look_forward_windows = 5
    window_count_threshold = 4

    # Always trim off the first few values as these are often dodgy.
    pos = initial_trim_size

    # Look at the stdev of the signal in the upcoming windows. Trimming is finished when:
    #  1. the next window has a high stdev
    #  2. enough of the other upcoming windows have a high stdev
    while True:
        next_window_stdev = get_window_stdev(signal, pos, 0, trim_increment)
        if next_window_stdev > stdev_threshold:
            upcoming_window_stdevs = [get_window_stdev(signal, pos, i, trim_increment)
                                      for i in range(look_forward_windows)]
            num_high_stdevs = sum(1 if x > stdev_threshold else 0
                                  for x in upcoming_window_stdevs)
            if num_high_stdevs >= window_count_threshold:
                return pos
        pos += trim_increment


class CannotTrim(IndexError):
    pass


def get_window_stdev(signal, current_pos, window_num, increment):
    window_start = current_pos + (window_num * increment)
    window_end = window_start + increment
    if window_end > len(signal):
        raise CannotTrim
    return np.std(signal[window_start:window_end])


if __name__ == '__main__':
    main()
