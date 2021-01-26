#!/usr/bin/env bash

# Set DPI
dpi=300

mkdir -p PNG_figures/

echo "Plotting Figure 9, panel B..."
python plot-staircase.py --dpi "$dpi"

echo "Plotting Figure 9, panel C..."
python plot-staircase.py --mutant 2 --dpi "$dpi"

echo "Plotting Figure 9, panel D..."
python plot-params-control.py --dpi "$dpi"
python plot-params.py --dpi "$dpi"
python plot-params.py --mutant 2 --dpi "$dpi"

mkdir -p PNG_figures/Figure9/
mv PNG_figures/*.png PNG_figures/Figure9/