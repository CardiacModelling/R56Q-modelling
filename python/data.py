#!/usr/bin/env python3
#
# Python module that knows where all the data, models, and protocols are, and
# can load them.
#
from __future__ import division, print_function
import inspect
import myokit
import numpy as np
import os

# Get root of this project
try:
    frame = inspect.currentframe()
    ROOT = os.path.dirname(inspect.getfile(frame))
finally:
    del(frame) # Must be manually deleted
ROOT = os.path.join(ROOT, '..')


# Data directory
DATA = os.path.join(ROOT, 'data')

# Model directory
MODEL = os.path.join(ROOT, 'model-and-protocols')

# Protocol directory
PROTO = os.path.join(ROOT, 'model-and-protocols')


def load(cell, protocol, mutant_str, cached=None):
    """
    Returns data for the given cell and protocol, with capacitance filtering
    applied.

    Arguments:

    ``cell``
        The cell to use (integer).
    ``protocol``
        The protocol to use (integer)
    ``cached``
        Optional cached data. If given, this will be returned directly.

    Returns a myokit DataLog.
    """
    if cached is not None:
        return cached

    # Get path to data file
    data_files = {
        1: os.path.join(DATA, 'SW/sine-wave-' + mutant_str + '-cell-' + str(cell)),
        2: os.path.join(DATA, 'Staircase/staircase1-' + mutant_str + '-cell-' + str(cell)),
        3: os.path.join(DATA, 'Staircase/staircase2-' + mutant_str + '-cell-' + str(cell)),
        4: os.path.join(DATA, 'Activation/activation-' + mutant_str + '-cell-' + str(cell)),
        5: os.path.join(DATA, 'Inactivation/inactivation-' + mutant_str + '-cell-' + str(cell)),
        6: os.path.join(DATA, 'AP/complex-AP-' + mutant_str + '-cell-' + str(cell)),
        7: os.path.join(DATA, 'Conductance/conduct-' + mutant_str + '-cell-' + str(cell))
    }
    data_file = data_files[protocol]

    # Load protocol for capacitance filtering.
    protocol = load_myokit_protocol(protocol)

    # Load data from zip or csv
    if os.path.exists(data_file + '.zip'):
        print('Loading ' + data_file + '.zip')
        log = myokit.DataLog.load(data_file + '.zip').npview()
    else:
        print('Loading ' + data_file + '.csv')
        log = myokit.DataLog.load_csv(data_file + '.csv').npview()
        log.save(data_file + '.zip')

    # Apply capacitance filtering
    dt = 0.1
    signals = [log.time(), log['current']]
    voltage = 'voltage' in log
    if voltage:
        signals.append(log['voltage'])
    signals = capacitance(protocol, dt, *signals)

    log = myokit.DataLog()
    log.set_time_key('time')
    log['time'] = signals[0]
    log['current'] = signals[1]
    if voltage:
        log['voltage'] = signals[2]

    # Return
    return log


def load_model(which_model):
    """
    Loads the selected model.
    """
    if which_model == 1:
        return myokit.load_model(os.path.join(MODEL, 'CCOI-ikr-markov-voffset.mmt'))
    else:
        return myokit.load_model(os.path.join(MODEL, 'M10-ikr-markov-voffset.mmt'))


def load_myokit_protocol(protocol):
    """
    Loads the Myokit protocol with the given index (1-8).
    """
    protocol_files = {
        1: os.path.join(PROTO, 'sine-wave-ramp-steps.mmt'),
        2: os.path.join(PROTO, 'staircase1.mmt'),
        3: os.path.join(PROTO, 'staircase2.mmt'),
        4: os.path.join(PROTO, 'activation-ramp.mmt'),
        5: os.path.join(PROTO, 'inactivation-ramp.mmt'),
        6: os.path.join(PROTO, 'pr6-ap-steps.mmt'),
        7: os.path.join(PROTO, 'conductance.mmt'),
        8: os.path.join(PROTO, 'staircase-within-staircase.mmt')
    }

    protocol = protocol_files[protocol]

    # Load Myokit protocol
    return myokit.load_protocol(protocol)


def load_ap_protocol():
    """
    Returns a tuple ``(times, values)`` representing AP protocol.
    """
    data_file = os.path.join(DATA, 'Protocols', 'AP-protocol')

    # Load data from zip or csv
    if os.path.exists(data_file + '.zip'):
        print('Loading ' + data_file + '.zip')
        log = myokit.DataLog.load(data_file + '.zip').npview()
    else:
        print('Loading ' + data_file + '.csv')
        log = myokit.DataLog.load_csv(data_file + '.csv').npview()
        log.save(data_file + '.zip')

    return log


def load_RPR_protocol(cell, mutant_str):
    """
    Returns a tuple ``(times, values)`` representing PPR protocol.
    """
    if cell in {4, 5} and mutant_str in {'WT', 'WT-RPR'}:
        data_file = os.path.join(DATA, 'Protocols', 'RPR-protocol-alt')
    else:
        data_file = os.path.join(DATA, 'Protocols', 'RPR-protocol')

    # Load data from zip or csv
    if os.path.exists(data_file + '.zip'):
        print('Loading ' + data_file + '.zip')
        log = myokit.DataLog.load(data_file + '.zip').npview()
    else:
        print('Loading ' + data_file + '.csv')
        log = myokit.DataLog.load_csv(data_file + '.csv').npview()
        log.save(data_file + '.zip')

    return log


def load_protocol_values(protocol, mutant_str, cell, no_filter=True, which_model=1):
    """
    Returns a (capacitance filtered) tuple ``(times, voltages)`` for the
    selected ``protocol``.
    """
    p = load_myokit_protocol(protocol)

    if protocol == 1:
        m = load_model(which_model)
        m.get('membrane.V').set_rhs(
            'piecewise(engine.time >= 300.0 and engine.time < 900.0,'
            + ' -140'
            + ' + 0.1 * (engine.time - 300.0),'
            + 'engine.time >= 3599.9 and engine.time < 7100.1,'
            + ' - 30'
            + ' + 54 * sin(0.007 * (engine.time - 3100.1))'
            + ' + 26 * sin(0.037 * (engine.time - 3100.1))'
            + ' + 10 * sin(0.190 * (engine.time - 3100.1)),'
            + 'engine.time >= 7100.1 and engine.time < 7200.1,'
            + ' -70'
            + ' - 0.3 * (engine.time - 7100.0)'
            + ', engine.pace)')
        p = load_myokit_protocol(protocol)
    elif protocol == 2:
        m = load_model(which_model)
        if cell in {4, 5} and mutant_str in {'WT', 'WT-RPR'}:
            m.get('membrane.V').set_rhs(
                'piecewise(engine.time <= 1296.8,'
                + ' 0,'
                + 'engine.time >= 14410.1 and engine.time < 14510.0,'
                + ' -70'
                + ' - 0.4 * (engine.time - 14410.1)'
                + ', engine.pace)')
        else:
            m.get('membrane.V').set_rhs(
               'piecewise(engine.time <= 1236.2,'
               + ' -80,'
               + 'engine.time >= 14410.1 and engine.time < 14510.0,'
               + ' -70'
               + ' - 0.4 * (engine.time - 14410.1)'
               + ', engine.pace)')

        p = load_myokit_protocol(protocol)
    elif protocol == 3:
        m = load_model(which_model)
        if cell in {4, 5} and mutant_str in {'WT', 'WT-RPR'}:
            m.get('membrane.V').set_rhs(
                'piecewise(engine.time > 300.1 and engine.time <= 700.1,'
                + ' -120'
                + ' + 0.1 * (engine.time - 300.1),'
                + 'engine.time >= 14410.1 and engine.time < 14478.8,'
                + ' -70'
                + ' - 0.4 * (engine.time - 14410.1),'
                + 'engine.time >= 14478.8,'
                + ' 0'
                + ', engine.pace)')
        else:
            m.get('membrane.V').set_rhs(
                'piecewise(engine.time > 300.1 and engine.time <= 700.1,'
                + ' -120'
                + ' + 0.1 * (engine.time - 300.1),'
                + 'engine.time >= 13000.0,'
                + ' -80,'
                + 'engine.pace)')
        p = load_myokit_protocol(protocol)
        s = myokit.Simulation(m, p)
        tmax = p.characteristic_time()
        t = np.arange(0, tmax, 0.1)
        v = s.run(tmax, log=['membrane.V'], log_times=t)
        v = np.array(v['membrane.V'])
    elif protocol in {4,5}:
        m = load_model(which_model)
        t = np.arange(0, p.characteristic_time(), 0.1)
        v = np.array(p.value_at_times(t))
        s = myokit.Simulation(m, p)
    elif protocol == 6:
        m = load_model(which_model)
        log = load_ap_protocol().npview()
        t, v = log['time'], log['voltage']
        s = myokit.Simulation(m, p)
        s.set_fixed_form_protocol(t, v)
    elif protocol == 8:
        m = load_model(which_model)
        m.get('membrane.V').set_rhs(
           'piecewise(engine.time <= 800.0,'
           + ' -80,'
           + 'engine.time >= 14410.1 and engine.time < 14510.0,'
           + ' -70'
           + ' - 0.4 * (engine.time - 14410.1)'
           + ', engine.pace)')
        p = load_myokit_protocol(protocol)
        s = myokit.Simulation(m, p)
        tmax = p.characteristic_time()
        t = np.arange(0, tmax, 0.1)
        v = s.run(tmax, log=['membrane.V'], log_times=t)
        v = np.array(v['membrane.V'])    
    else:
        t = np.arange(0, p.characteristic_time(), 0.1)
        v = np.array(p.value_at_times(t))

    if no_filter:
        return p
    else:
        return capacitance(p, 0.1, t, v)


def capacitance(protocol, dt, *signals):
    """
    Creates and applies a capacitance filter, based on a Myokit protocol.

    Arguments:

    ``protocol``
        A Myokit protocol.
    ``dt``
        The sampling interval of the given signals.
    ``signals``
        One or more signal files to filter.

    Returns a filtered version of the given signals.
    """
    cap_duration = 1    # Kylie and Michael used 5 ms, this is far too long for 37C data, so 1 ms is used instead
    fcap = np.ones(len(signals[0]), dtype=int) # fcap = length of signals[0], which is time
    steps = [step for step in protocol]
    for step in steps[1:]:
        i1 = int(step.start() / dt)
        i2 = i1 + int(cap_duration / dt)
        fcap[i1:i2] = 0
    fcap = fcap > 0

    debug = False
    if debug:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(signals[0], signals[1])
        plt.plot(signals[0], signals[1])
        for step in steps[1:]:
            plt.axvline(step.start(), color='green', alpha=0.25)
        plt.show()

    # Apply filter
    return [x[fcap] for x in signals]


def model_path(model_file):
    """
    Returns the path to the given Myokit model file.
    """
    return os.path.join(MODEL, model_file)


def protocol_path(protocol_file):
    """
    Returns the path to the given Myokit protocol file.
    """
    return os.path.join(PROTO, protocol_file)

