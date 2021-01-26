#!/usr/bin/env bash

# Set DPI
dpi=300

mkdir -p PNG_figures/

echo "Plotting Figure 13, panel A..."
python plot-ORd.py --dpi "$dpi"

echo "Plotting Figure 13, panel B..."
python plot-ORd.py --RPR --dpi "$dpi"

echo "Plotting Figure 13, panel C..."
python plot-ORd-S1S2.py --dpi "$dpi"

mkdir -p PNG_figures/Figure13/
mv PNG_figures/*.png PNG_figures/Figure13/