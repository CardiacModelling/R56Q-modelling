#!/usr/bin/env bash

# Settings
dpi=300
fig_num=9

mkdir -p PNG_figures/

echo "Plotting Figure $fig_num, panel B..."
python plot-staircase.py --dpi "$dpi"

echo "Plotting Figure $fig_num, panel C..."
python plot-staircase.py --mutant 2 --dpi "$dpi"

echo "Plotting Figure $fig_num, panel D..."
python plot-params-control.py --dpi "$dpi"
python plot-params.py --dpi "$dpi"
python plot-params.py --mutant 2 --dpi "$dpi"

mkdir -p PNG_figures/Figure"$fig_num"/
mv PNG_figures/*.png PNG_figures/Figure"$fig_num"/
