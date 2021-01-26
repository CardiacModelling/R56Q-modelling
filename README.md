# Modelling the effects of the R56Q hERG channel mutation

This repository contains scripts for optimising models of the cardiac hERG ion channel to wild type (WT) and R56Q mutant channel electrophysiology data using [Myokit](http://myokit.org) and [PINTS](https://github.com/pints-team/pints) modules in Python.
This code is associated with the paper:

_"The hERG channel activator, RPR260243, reduces arrhythmogenicity of the R56Q LQTS mutation by enhancing repolarizing drive in the refractory period_". (Currently under review). Kemp, J. M., Whittaker, D. G., Venkateshappa, R., Pang, Z., Johal, R., Sergeev, V., Tibbits, G. F., Mirams, G. R., Claydon, T. W.

## Installation

It is recommended to install libraries and run scripts in a virtual environment to avoid version conflicts between different projects. In order to do this, follow these steps:
- git clone the repository
- Set up a virtual environment using `virtualenv venv` or if you have both python 2 and 3: `virtualenv --python=python3 venv`. If that doesn't work you may need to install virtualenv first: `pip install virtualenv`.
- Activate the virtual environment using `source venv/bin/activate`. Simply type `deactivate` to exit the virtual environment at the end of a session.
- Install the required packages by typing `pip install -r requirements.txt`

## Protocols and data

### Voltage protocols

For calibration and validation of the mathematical model, this study used a novel collection of protocols (shown below), including the previously-described 'staircase' protocol by [(Lei _et al._, 2019)](https://www.sciencedirect.com/science/article/pii/S0006349519305971).
<img src="https://github.com/CardiacModelling/R56Q-modelling/blob/main/figures/Paper_figures/full-protocol.png" height="400">

### Experimental data

## Running

### Parameter inference

In order to reproduce the parameter fitting, simply run `python cmaesfit.py` for WT, and `python cmaesfit.p --mutant 2` for R56Q. This should take several hours running in parallel, so is best performed using a HPC resource.

### Plotting results

## Acknowledging this work

If you publish any work based on the contents of this repository please cite:

[PLACEHOLDER]