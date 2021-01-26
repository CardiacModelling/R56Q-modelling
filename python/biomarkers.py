#!/usr/bin/env python3
#
# Biomarker based functions
#
from __future__ import division
import numpy as np

def steady_state_activation(pr_steps, pr_voltages, log=None, normalise=False):
    """
    Returns the steady state of activation.

    Arguments:

    ``log``
        An optional datalog with the data for the given cell.

    Returns a tuple ``(voltages, activation)`` where ``voltages`` and
    ``activation`` are lists of equal lengths.
    """
    # Load data, or use cached
    time = log['time']
    current = log['current']

    steps = pr_steps
    voltages = pr_voltages

    # Find peaks
    cpeaks = []
    for i, j in steps:
        c = current[i:i + j]
        z = np.argmax(c)
        zlo = max(z - 5, 0)
        zhi = zlo + 10
        cpeaks.append(np.mean(c[zlo:zhi]))

    # Normalise
    # Note: Division by (V - E) not needed; all measured at same voltage!
    if normalise:
        cpeaks /= np.max(cpeaks)
    else:
        cpeaks = np.array(cpeaks)

    return voltages, cpeaks

def iv_curve_activation(pr_steps, pr_voltages, log=None, normalise=False):
    """
    Returns the steady state of activation.

    Arguments:

    ``log``
        An optional datalog with the data for the given cell.

    Returns a tuple ``(voltages, activation)`` where ``voltages`` and
    ``activation`` are lists of equal lengths.
    """
    # Load data, or use cached
    time = log['time']
    current = log['current']

    steps = pr_steps
    voltages = pr_voltages

    # Find peaks
    cpeaks = []
    for i, j in steps:
        c = current[i-j:i]
        z = np.argmax(c)
        zlo = max(z - 5, 0)
        zhi = zlo + 10
        cpeaks.append(np.mean(c[zlo:zhi]))

    # Normalise
    # Note: Division by (V - E) not needed; all measured at same voltage!
    if normalise:
        cpeaks /= np.max(cpeaks)
    else:
        cpeaks = np.array(cpeaks)

    return voltages, cpeaks

def steady_state_inactivation(
        pr_steps, pr_voltages, log=None, include_minus_90=False, estimate_erev=False, erev=-88, normalise=False):
    """
    Returns the steady state of inactivation and an IV curve

    Arguments:

    ``log``
        An optional datalog with the data for the given cell.

    Returns a tuple ``(v1, inactivation, v2, iv_curve)`` where ``v1`` and
    ``inactivation`` are of equal length, and ``v2`` and ``iv_curve`` are of
    equal length.
    """
    # Load data, or use cached
    time = log['time']
    current = log['current']

    steps = pr_steps
    voltages = np.array(pr_voltages)

    # Find peaks
    peaks = []
    for k, step in enumerate(steps):
        i, j = step
        c = current[i:i + j]
        peaks.append(current[i + np.argmax(np.abs(c))])
    peaks = np.array(peaks)

    if estimate_erev:
        # Estimate reversal potential from IV curve
        irev = np.argmax(peaks >= 0)
        x1, x2 = voltages[irev - 1], voltages[irev]
        y1, y2 = peaks[irev - 1], peaks[irev]
        erev = x1 - y1 * (x2 - x1) / (y2 - y1)

    # Calculate steady state
    v1 = voltages
    g = peaks / (voltages - erev)
    if not include_minus_90:
        # Remove point at -90, too close to Erev!
        v1 = np.concatenate((voltages[:2], voltages[3:]))
        g = np.concatenate((g[:2], g[3:]))
    if normalise:
        g /= np.max(g)

    # Return
    return v1, g, voltages, peaks

def time_constant_of_inactivation(pr_steps, pr_voltages, log=None):
    """
    Returns time constants of inactivation.

    Returns a tuple ``(voltages, time_constants)`` where ``voltages`` and
    ``time_constants`` are lists of equal lengths.
    """
    # Load data, or use cached
    time = log['time']
    current = log['current']

    steps = pr_steps
    voltages = pr_voltages

    debug = False
    if debug:
        import matplotlib.pyplot as plt
        plt.figure()

    # Curve fitting, single exponential
    from scipy.optimize import curve_fit

    def f(t, a, b, c):
        if c <= 0:
            return np.ones(t.shape) * float('inf')
        return a + b * np.exp(-t / c)

    taus = []
    for k, step in enumerate(steps):
        i, j = step
        imin = i + np.argmin(current[i:i + j])
        t = time[i:imin] - time[i]
        c = current[i:imin]

        # Guess some parameters and fit
        p0 = current[imin - 1], current[i] - current[imin - 1], 10
        popt, pcov = curve_fit(f, t, c, p0)
        taus.append(popt[2])

        if debug:
            print(voltages[k], popt[2])
            plt.plot(t, c)
            plt.plot(t, f(t, *popt), 'k:')

    if debug:
        plt.show()

    return voltages, taus

def time_constant_of_deactivation(pr_steps, pr_voltages, log=None, erev=-88):
    """
    Returns time constants of inactivation.

    Arguments:

    ``pr4_log``
        An optional datalog with the data for the given cell.

    Returns a tuple ``(voltages, time_constants)`` where ``voltages`` and
    ``time_constants`` are lists of equal lengths.
    """
    # Load data, or use cached
    time = log['time']
    current = log['current']

    steps = pr_steps
    voltages = pr_voltages

    debug = False
    if debug:
        import matplotlib.pyplot as plt
        plt.figure()

    # Curve fitting, single exponential
    from scipy.optimize import curve_fit

    def f(t, a, b, c):
        if c <= 0:
            return np.ones(t.shape) * float('inf')
        return a + b * np.exp(-t / c)

    taus = []
    for k, step in enumerate(steps):
        i, j = step
        v = voltages[k]
        if v < erev:
            ipeak = i + np.argmin(current[i:i + j])
        else:
            ipeak = i + np.argmax(current[i:i + j])
        t = time[ipeak:ipeak+1000] - time[ipeak]
        c = current[ipeak:ipeak+1000]

        # Guess some parameters and fit
        p0 = current[j - 1], current[ipeak] - current[j - 1], 10
        popt, pcov = curve_fit(f, t, c, p0, maxfev=8000)
        taus.append(popt[2])

        if debug:
            print(voltages[k], popt[2])
            plt.plot(t, c)
            plt.plot(t, f(t, *popt), 'k:')

    if debug:
        plt.show()

    return voltages, taus
