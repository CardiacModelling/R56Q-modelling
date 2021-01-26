#!/usr/bin/env bash

# Settings
dpi=300
fig_num=11

mkdir -p PNG_figures/

echo "Plotting Figure $fig_num, panel A..."
python plot-inactivation.py --dpi "$dpi"

echo "Plotting Figure $fig_num, panel B..."
python plot-inactivation.py --mutant 2 --dpi "$dpi"

mkdir -p PNG_figures/Figure"$fig_num"/
mv PNG_figures/*.png PNG_figures/Figure"$fig_num"/
