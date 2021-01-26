# Modelling the effects of the R56Q hERG channel mutation

This repository contains scripts for optimising models of the cardiac hERG ion channel to wild type (WT) and R56Q mutant hERG channel electrophysiology data using [Myokit](http://myokit.org) and [PINTS](https://github.com/pints-team/pints) modules in Python.
This code is associated with the paper:

_"The hERG channel activator, RPR260243, reduces arrhythmogenicity of the R56Q LQTS mutation by enhancing repolarizing drive in the refractory period_" (Currently under review). Kemp, J. M., Whittaker, D. G., Venkateshappa, R., Pang, Z., Johal, R., Sergeev, V., Tibbits, G. F., Mirams, G. R., Claydon, T. W.

## Installation

The codes in this repository use Python 3. It is recommended to install libraries and run scripts in a virtual environment to avoid version conflicts between different projects. In order to do this, follow these steps:
- git clone the repository
- Set up a virtual environment using `virtualenv venv` or if you have both Python 2 and 3: `virtualenv --python=python3 venv`. If that doesn't work you may need to install virtualenv first: `pip install virtualenv`
- Activate the virtual environment using `source venv/bin/activate`. Simply type `deactivate` to exit the virtual environment at the end of a session.
- Install the required packages by typing `pip install -r requirements.txt`

## Protocols and data

### Voltage protocols

For calibration and validation of the mathematical model, the study used a novel collection of protocols (shown below), including the previously-described 'staircase' protocol by [(Lei _et al._, 2019)](https://www.sciencedirect.com/science/article/pii/S0006349519305971).

<img src="https://github.com/CardiacModelling/R56Q-modelling/blob/main/figures/Paper_figures/full-protocol.png">

The voltage clamp waveform itself can be found in `.csv` form in [data/Protocols/](https://github.com/CardiacModelling/R56Q-modelling/tree/main/data/Protocols) as `RPR-protocol.csv`, which contains a list of voltages (in mV) and times (in ms) that comprise the protocol at a sampling frequency of 10 kHz.

### Experimental data

Experimental data obtained from the above electrophysiology protocol (recorded at 37C in HEK cells) have been processed, separated into constituent protocols, and stored in folders in [data/](https://github.com/CardiacModelling/R56Q-modelling/tree/main/data). The data are named according to `[protocol]-[mutant]-cell-[cell number].csv`, and are stored in sub-folders according to protocol. The activation protocol data for WT, cell 2 can thus be found in `activation-WT-cell-2.csv` in the [data/Activation/](https://github.com/CardiacModelling/R56Q-modelling/blob/main/data/Activation) folder.

## Running

### Parameter inference

The parameter inference in this study uses the [CMA-ES](https://www.mitpressjournals.org/doi/abs/10.1162/106365603321828970) algorithm to search the parameter space. In order to reproduce the parameter fitting, simply run
- `python cmaesfit.py` and
- `python cmaesfit.p --mutant 2`

for WT and R56Q, respectively. Each of these could take several hours running in parallel, so are best performed using a HPC resource. The number of repeats can be reduced from the default 25 using the `--repeats` input argument, although this may not explore the parameter space sufficiently to find a global minimum. The default model used, described as `M10` throughout the code, corresponds to the structure shown below (also in manuscript Figure 9A):

<img src="https://github.com/CardiacModelling/R56Q-modelling/blob/main/figures/Paper_figures/markov-chain.png">

Once finished, the results of the CMA-ES optimisation are deposited in [cmaesfits/](https://github.com/CardiacModelling/R56Q-modelling/tree/main/cmaesfits). The maximum likelihood estimate parameters for the default runs can be found in `WT-model-M10-fit-staircase1-artefact.txt` and `R56Q-model-M10-fit-staircase1-artefact.txt` for WT and R56Q, respectively. Detailed breakdowns of the scores and parameters for each repeat can be found in the accompanying `parameter` and `log` files.

### Plotting results

Most of the plots found in the manuscript are composite plots consisting of multiple, individually-generated panels. Nonetheless, these can be generated easily using python scripts provided in [figures/](https://github.com/CardiacModelling/R56Q-modelling/tree/main/figures). It is recommended to run the scripts with `--show` to simply check the results on screen, _or_ with `--dpi 300` in order to generate high-quality, publication-standard figures in `.png` format (stored in sub-folders of the [figures/](https://github.com/CardiacModelling/R56Q-modelling/tree/main/figures) directory). Where results from a single cell are shown, these correspond to cell 2 for the WT recordings, and cell 5 for R56Q (which can be identified based on the filename). In order to quickly, generate all panels for a given figure, bash scripts are provided. If these cannot be executed, simply type `chmod u+x [bash script name].sh` or, alternatively, `chmod u+x *.sh`.

- Manuscript Figure 9A is a schematic drawing. To generate Figure panels 9B-D, run `./figure9-panels.sh`
- To generate panels for Figures 10-13, run `./figure[figure number]-panels.sh`

All final, composite figure files (pertaining to the modelling) used in the manuscript are stored in [figures/Paper_figures/](https://github.com/CardiacModelling/R56Q-modelling/tree/main/figures/Paper_figures).

## Acknowledging this work

If you publish any work based on the contents of this repository please cite:

[PLACEHOLDER]