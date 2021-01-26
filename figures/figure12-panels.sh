#!/usr/bin/env bash

# Set DPI
dpi=300

mkdir -p PNG_figures/

echo "Plotting Figure 12, panel A..."
python plot-AP.py --dpi "$dpi"

echo "Plotting Figure 12, panel B..."
python plot-AP.py --mutant 2 --dpi "$dpi"

echo "Plotting Figure 12, panel C..."
python plot-AP-compare.py --dpi "$dpi"

mkdir -p PNG_figures/Figure12/
mv PNG_figures/*.png PNG_figures/Figure12/