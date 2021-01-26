#!/usr/bin/env bash

# Settings
dpi=300
fig_num=12

mkdir -p PNG_figures/

echo "Plotting Figure $fig_num, panel A..."
python plot-AP.py --dpi "$dpi"

echo "Plotting Figure $fig_num, panel B..."
python plot-AP.py --mutant 2 --dpi "$dpi"

echo "Plotting Figure $fig_num, panel C..."
python plot-AP-compare.py --dpi "$dpi"

mkdir -p PNG_figures/Figure"$fig_num"/
mv PNG_figures/*.png PNG_figures/Figure"$fig_num"/
