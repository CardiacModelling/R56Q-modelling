#!/usr/bin/env bash

# Set DPI
dpi=300

mkdir -p PNG_figures/

echo "Plotting Figure 11, panel A..."
python plot-inactivation.py --dpi "$dpi"

echo "Plotting Figure 11, panel B..."
python plot-inactivation.py --mutant 2 --dpi "$dpi"

mkdir -p PNG_figures/Figure11/
mv PNG_figures/*.png PNG_figures/Figure11/