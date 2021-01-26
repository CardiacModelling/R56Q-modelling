#!/usr/bin/env bash

# Settings
dpi=300
fig_num=13

mkdir -p PNG_figures/

echo "Plotting Figure $fig_num, panel A..."
python plot-ORd.py --dpi "$dpi"

echo "Plotting Figure $fig_num, panel B..."
python plot-ORd.py --RPR --dpi "$dpi"

echo "Plotting Figure $fig_num, panel C..."
python plot-ORd-S1S2.py --dpi "$dpi"

mkdir -p PNG_figures/Figure"$fig_num"/
mv PNG_figures/*.png PNG_figures/Figure"$fig_num"/
