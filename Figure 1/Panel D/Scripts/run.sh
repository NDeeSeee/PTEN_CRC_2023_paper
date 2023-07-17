#!/bin/bash

PYTHON_PATH="python"
SCRIPT_DIR="scripts"

 
mkdir -p results

# find linear hotspots
# arguments order: mutations frequency | p-value cutoff
$PYTHON_PATH $SCRIPT_DIR/generate_linear_hotspots.py data/IG_muts.tsv 0.005 > results/linear_hotspots.tsv

# find 3D hotspots
# arguments order: 3d interactions | mutations frequency | table with random indicies | linear hotspots | output file
$PYTHON_PATH $SCRIPT_DIR/generate_3d_hotspots.py data/PTEN_interactions.tsv data/IG_muts.tsv data/randomTable.tsv results/linear_hotspots.tsv results/3dHotspots.tsv

# Compile to one file
$PYTHON_PATH $SCRIPT_DIR/prepare_results.py data/IG_muts.tsv results/3dHotspots.tsv results.xlsx





# # find sliding window enrichments
# # arguments order: mutations frequency | p-value cutoff
# $PYTHON_PATH $SCRIPT_DIR/generate_sliding_windows.py data/variants_2019_nras_counts.tsv results/sliding_window.tsv


