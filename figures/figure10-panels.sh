#!/usr/bin/env bash

# Set DPI
dpi=300

mkdir -p PNG_figures/

echo "Plotting Figure 10, panel A..."
python plot-activation.py --dpi "$dpi"

echo "Plotting Figure 10, panel B..."
python plot-activation.py --mutant 2 --dpi "$dpi"

mkdir -p PNG_figures/Figure10/
mv PNG_figures/*.png PNG_figures/Figure10/