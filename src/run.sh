#!/bin/sh
# Created: Dec, 17, 2025 16:27:21 by Wataru Fukuda
set -eu

# BASE=$(readlink -f $(dirname $0))
# N=(10 20 30 50 80 120 160 210 270 350)
N=(10 20 30 50 80 120 160 210 270 350 450 600 800 1000)
# N=(350)

matlab=/Applications/MATLAB_R2024b.app/bin/matlab

for i in "${N[@]}"; do
  export N_INPUT=$i
  $matlab -batch WBM_Bspline
done

