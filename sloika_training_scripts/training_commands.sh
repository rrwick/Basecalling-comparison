#!/usr/bin/env bash

# Copyright 2019 Ryan Wick (rrwick@gmail.com)
# https://github.com/rrwick/Basecalling-comparison

# This program is free software: you can redistribute it and/or modify it under the terms of the
# GNU General Public License as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version. This program is distributed in the hope that it
# will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You
# should have received a copy of the GNU General Public License along with this program. If not,
# see <http://www.gnu.org/licenses/>.


# Edit the following paths before running, as appropriate for your environment:
out_dir=training_output
training_input_dir=/path/to/training_data       # this directory must contain the hdf5 files made by Sloika chunkify
sloika_dir=/path/to/sloika_dir                  # directory of Sloika clone (https://github.com/rrwick/sloika)
model=/path/to/sloika/models/raw_0.98_rgrgr.py  # model file - change this to change the neural network architecture


# Do the actual training:
"$sloika_dir"/bin/train_network.py raw --min_prob 1e-5 --reload_after_batches 500 --save_every 500 --niteration 100000 "$model" "$out_dir" "$training_input_dir"/*.hdf5

# Convert the trained models into both CPU and GPU appropriate forms:
for f in "$out_dir"/model_*.pkl; do
    "$sloika_dir"/misc/model_convert.py --target cpu $f "$f".cpu
    "$sloika_dir"/misc/model_convert.py --target gpu $f "$f".gpu
done
