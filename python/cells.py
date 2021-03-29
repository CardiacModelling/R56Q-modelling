#!/usr/bin/env python3
#
from __future__ import division, print_function
import numpy as np


def lower_conductance(cell, mutant_str):
    """
    Returns a lower limit for the conductance of the cell with the
    given integer index ``cell``.
    """
    #
    # Guesses for lower conductance
    #
    if mutant_str == 'WT':
        lower_conductances = {
            1: 0.114,
            2: 0.108,
            3: 0.151,
            4: 0.169,
            5: 0.0821
        }
    elif mutant_str == 'R56Q':
        lower_conductances = {
            1: 0.0938,
            2: 0.125,
            3: 0.0732,
            4: 0.0642,
            5: 0.154,
            6: 0.0945
        }
    elif mutant_str == 'WT-RPR':
        lower_conductances = {
            1: 0.0590,
            2: 0.0876,
            3: 0.110,
            4: 0.0989,
            5: 0.0422
        }
    else:
        lower_conductances = {
            1: 0.0522,
            2: 0.0219,
            3: 0.0520,
            4: 0.0302,
            5: 0.0623,
            6: 0.0444
        }
    return lower_conductances[cell]


def reversal_potential(temperature):
    """
    Calculates the reversal potential for Potassium ions, using the Nernst
    equation for a given ``temperature`` in degrees celsius and the internal
    and external [K]+ concentrations used in the experiments.
    """
    T = 273.15 + temperature
    F = 96485
    R = 8314
    K_i = 130
    k_o = 4
    return ((R*T)/F) * np.log(k_o/K_i)


def temperature(cell):
    """
    Returns the temperature (in degrees Celsius) for the given integer index
    ``cell``.
    """
    temperatures = {
        0: 21.0,
    }
    return temperatures[cell]

def ek_computed():
    return -93.04

def ek_WT(cell):
    """
    Returns the WT reversal potential (in mV) for the given integer index
    ``cell``.
    """
    reversal_potentials = {
        1: -91.6,
        2: -92.8,
        3: -95.1,
        4: -92.3,
        5: -106.1
    }
    return reversal_potentials[cell]


def ek_WT_RPR(cell):
    """
    Returns the WT-RPR reversal potential (in mV) for the given integer index
    ``cell``.
    """
    reversal_potentials = {
        1: -87.4,
        2: -92.1,
        3: -96.1,
        4: -93.1,
        5: -106.1
    }
    return reversal_potentials[cell]


def ek_R56Q(cell):
    """
    Returns the R56Q reversal potential (in mV) for the given integer index
    ``cell``.
    """
    reversal_potentials = {
        1: -96.0,
        2: -95.0,
        3: -90.5,
        4: -94.5,
        5: -94.5,
        6: -101.0
    }
    return reversal_potentials[cell]

def ek_R56Q_RPR(cell):
    """
    Returns the R56Q-RPR reversal potential (in mV) for the given integer index
    ``cell``.
    """
    reversal_potentials = {
        1: -92.2,
        2: -89.0,
        3: -90.1,
        4: -95.0,
        5: -93.1,
        6: -93.4
    }
    return reversal_potentials[cell]


