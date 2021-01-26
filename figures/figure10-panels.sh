#!/usr/bin/env bash

# Settings
dpi=300
fig_num=10

mkdir -p PNG_figures/

echo "Plotting Figure $fig_num, panel A..."
python plot-activation.py --dpi "$dpi"

echo "Plotting Figure $fig_num, panel B..."
python plot-activation.py --mutant 2 --dpi "$dpi"

mkdir -p PNG_figures/Figure"$fig_num"/
mv PNG_figures/*.png PNG_figures/Figure"$fig_num"/
