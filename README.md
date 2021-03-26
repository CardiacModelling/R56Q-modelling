# Modelling the effects of the R56Q hERG channel mutation

This repository contains scripts for optimising models of the cardiac hERG ion channel to wild type (WT) and R56Q mutant hERG channel electrophysiology data using [Myokit](http://myokit.org) and [PINTS](https://github.com/pints-team/pints) modules in Python. This code is associated with the paper:

***"Electrophysiological characterization of the hERG R56Q LQTS variant and targeted rescue by the activator RPR260243"*** (Currently under review). Kemp, J. M., Whittaker, D. G., Venkateshappa, R., Pang, Z., Johal, R., Sergeev, V., Tibbits, G. F., Mirams, G. R., Claydon, T. W.

## Prerequisites
The instructions in this section assume that you are using a Linux distribution with `python3`, `pip` and `virutalenv` installed. For example, if you are using _Ubuntu 20.04.1 LTS_, running the command ```sudo apt update``` followed by ```sudo apt install virtualenv python3``` will install these packages.

It is recommended to install libraries and run this repository's scripts in a virtual environment to avoid version conflicts between different projects.
To create such a virtual environment open a terminal and do the following:
- Clone the repository (for example, by using `git clone https://github.com/CardiacModelling/R56Q-modelling.git`)
- Type `cd R56Q-modelling` to navigate inside the repository
- Set up a virtual environment using `virtualenv --python=python3 venv`
- Activate the virtual environment using `source venv/bin/activate`
- Install the required packages by typing `pip install -r requirements.txt`.

When you are finished working with the repository, type `deactivate` to exit the virtual environment. The virtual environment can be used again with the `activate` command as shown above.

## Protocols and data

### Voltage protocols

For calibration and validation of the mathematical model, the study used a novel collection of protocols (shown below), including the previously-described 'staircase' protocol by [(Lei _et al._, 2019)](https://www.sciencedirect.com/science/article/pii/S0006349519305971).

<img src="https://github.com/CardiacModelling/R56Q-modelling/blob/main/figures/Paper_figures/full-protocol.png">

The full voltage clamp waveform itself can be found as a comma-separated values (`.csv`) file in [data/Protocols/](https://github.com/CardiacModelling/R56Q-modelling/tree/main/data/Protocols) called `RPR-protocol.csv`. This contains a list of voltages (in mV) and times (in ms) that comprise the protocol at a sampling frequency of 10 kHz.

### Experimental data

Experimental data obtained from the above electrophysiology protocol (recorded at 37°C in HEK cells) have been processed, separated into constituent protocols, and stored in [data/](https://github.com/CardiacModelling/R56Q-modelling/tree/main/data). The data files are named according to `[protocol]-[mutant]-cell-[cell number].csv`, and are stored in sub-folders according to protocol. The activation protocol data for WT, cell 2 can thus be found in `activation-WT-cell-2.csv` in the [data/Activation/](https://github.com/CardiacModelling/R56Q-modelling/blob/main/data/Activation) folder.

## Running

### Parameter inference

The parameter inference in this study uses the [CMA-ES](https://www.mitpressjournals.org/doi/abs/10.1162/106365603321828970) algorithm to search the parameter space. To reproduce the parameter fitting, generate fits for WT and R56Q by running
- `python cmaesfit.py` and
- `python cmaesfit.py --mutant 2`.

Each of these could take several hours running in parallel, so are best performed using a HPC resource. The number of repeats can be reduced from the default 20 using the `--repeats` input argument, although this may not explore the parameter space sufficiently to find a global minimum. The default model used, described as `M10` throughout the code, corresponds to the structure shown below (taken from manuscript Figure 9). The code also supports inference for a C-C-O-I structure hERG model (labelled `CCOI` in the code), which can be toggled by changing the `--model` input argument.

<img src="https://github.com/CardiacModelling/R56Q-modelling/blob/main/figures/Paper_figures/markov-chain.png">

Once finished, the results of the CMA-ES optimisation are deposited in [cmaesfits/](https://github.com/CardiacModelling/R56Q-modelling/tree/main/cmaesfits). The maximum likelihood estimate parameters for the default runs can be found in `WT-model-M10-fit-staircase1-artefact.txt` and `R56Q-model-M10-fit-staircase1-artefact.txt` for WT and R56Q, respectively. Detailed breakdowns of the scores and parameters for each repeat can be found in the accompanying `parameter` and `log` `.txt` files. For convenience, the results of the inference used for the actual modelling data in the manuscript, along with detailed parameter and score logs, can be found already in [cmaesfits/](https://github.com/CardiacModelling/R56Q-modelling/tree/main/cmaesfits).

### Plotting results

Most of the plots found in the manuscript are composite plots consisting of multiple, individually-generated panels. Nonetheless, these can be generated easily using Python scripts provided in [figures/](https://github.com/CardiacModelling/R56Q-modelling/tree/main/figures). It is recommended to run the scripts with `--show` to simply check the results on screen, _or_ with `--dpi 300` in order to generate high-quality, publication-standard figures in `.png` format (stored in sub-folders of the [figures/](https://github.com/CardiacModelling/R56Q-modelling/tree/main/figures) directory). Where results from a single cell are shown, these correspond to cell 2 for the WT recordings, and cell 5 for R56Q (which can be identified based on the filename). In order to quickly generate all panels for a given figure, Bash scripts are provided. To use them, open a terminal in the figures directory and do the following:

- Type `chmod u+x *.sh` to ensure that the scripts are executable
- Manuscript Figure 9A is a schematic drawing. To generate Figure panels 9B-D, run `./figure9-panels.sh`
- Similarly, to generate panels for Figures 10-13, run `./figure[figure number]-panels.sh`

All final, composite figure files (pertaining to the _in silico_ modelling) used in the manuscript are stored in [figures/Paper_figures/](https://github.com/CardiacModelling/R56Q-modelling/tree/main/figures/Paper_figures).

## Acknowledging this work

If you publish any work based on the contents of this repository please cite:

[PLACEHOLDER]

### Related publications

The scientific approach and code structure in this work owe a lot to the following previous, related publications from our lab which may also be of interest:

Beattie, K. A., Hill, A. P., Bardenet, R., Cui, Y., Vandenberg, J. I., Gavaghan, D. J., de Boer, T. P., Mirams, G. R.
(2018).
[Sinusoidal voltage protocols for rapid characterisation of ion channel kinetics](https://doi.org/10.1113/JP275733).
_The Journal of Physiology_, 596, 1813–1828.

Clerx, M., Beattie, K. A., Gavaghan, D. J., Mirams, G. R.
(2019).
[Four ways to fit an ion channel model](https://doi.org/10.1016/j.bpj.2019.08.001).
_Biophysical Journal_, 117, 2420-2437.

Lei, C. L., Clerx, M., Gavaghan, D. J., Polonchuk, L., Mirams, G. R., Wang, K.
(2019).
[Rapid characterisation of hERG channel kinetics I: using an automated high-throughput system](https://doi.org/10.1016/j.bpj.2019.07.029).
_Biophysical Journal_, 117, 2438-2454.

Lei, C. L., Clerx, M., Whittaker, D. G., Gavaghan D. J., de Boer, T. P., Mirams, G. R.
(2020).
[Accounting for variability in ion current recordings using a mathematical model of artefacts in voltage-clamp experiments](https://doi.org/10.1098/rsta.2019.0348).
_Philosophical Transactions of the Royal Society A_, 378: 20190348.
