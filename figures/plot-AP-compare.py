import myokit
import myokit.pacing as pacing
import numpy as np
import matplotlib
import matplotlib.pyplot as pl
import myokit.lib.markov as markov
import pints
import argparse
import os
import sys

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

# Load project modules
sys.path.append(os.path.abspath(os.path.join('../', 'python')))
import data

# Check input arguments
parser = argparse.ArgumentParser(
    description='Plot model and experimental data')
parser.add_argument('--model', type=int, default=2, metavar='N',
                    help='model number : 1 for CCOI, 2 for M10')
parser.add_argument('--mutant', type=int, default=1, metavar='N',
                    help='mutant number : 1 for WT, 2 for R56Q')
parser.add_argument('-p', '--protocol', type=int, default=2, metavar='N',
                    help='which protocol is used to fit the data: 1 for sine wave, 2 for staircase #1')
parser.add_argument("--show", action='store_true',
                    help="whether to show figures instead of saving them",
                    default=False)
parser.add_argument('--dpi', type=int, default=100, metavar='N',
                    help='what DPI to use for figures (suggested 300 for publication quality)')
args = parser.parse_args()

# Get model
p = myokit.load_protocol('../model-and-protocols/pr6-ap-steps.mmt')

current = 'ikr.IKr'

if args.protocol == 1:
    protocol_str = 'sine-wave'
elif args.protocol == 2:
    protocol_str = 'staircase1'
elif args.protocol == 3:
    protocol_str = 'staircase2'
elif args.protocol == 4:
    protocol_str = 'activation'
elif args.protocol == 5:
    protocol_str = 'inactivation'
elif args.protocol == 6:
    protocol_str = 'complex-AP'
else:
    protocol_str = 'staircase1-conductance'

# Run simulation
dt = 0.1

mutants = [1, 2, 3, 4]

fig, (a0, a1) = pl.subplots(2, 1, gridspec_kw={'height_ratios': [1, 3]}, figsize=(10, 4), dpi=args.dpi, constrained_layout=True)
ds = []

for m in mutants:
    mutant = m

    if mutant in {1, 3}:
        mutant_str = 'WT'
        model = args.model
    else:
        mutant_str = 'R56Q'
        model = args.model

    if args.model == 1:
        m = myokit.load_model('../model-and-protocols/CCOI-ikr-markov-voffset.mmt')
        model_str = 'CCOI'
        states = [
            'ikr.c1',
            'ikr.c2',
            'ikr.i',
            'ikr.o'
        ]
    elif args.model == 2:
        m = myokit.load_model('../model-and-protocols/M10-ikr-markov-voffset.mmt')
        model_str = 'M10'
        states = [
            'ikr.c1',
            'ikr.c2',
            'ikr.i',
            'ikr.ic1',
            'ikr.ic2',
            'ikr.o'
        ]
    else:
        print('Invalid model')
        sys.exit()

    n_params = int(m.get('misc.n_params').value())

    # Set steady state potential
    ss_V = -80

    x_found = np.loadtxt('../cmaesfits/' + mutant_str + '-model-' + model_str + '-fit-' + protocol_str + '-artefact.txt', unpack=True)

    parameters = []
    for i in range(n_params):
        parameters.append('ikr.p'+str(i+1))

    d = [
        'engine.time',
        'membrane.V',
        'ikr.IKr'
        ]

    # Run simulation
    m.get('ikr.p'+str(n_params+1)).set_rhs(0.1)

    print('Updating model to steady-state for ' + str(ss_V) + 'mV.')
    m.get('membrane.V').set_label('membrane_potential')

    mm = markov.LinearModel(m, states, parameters, current)

    if mutant in {1, 2}:
        x = mm.steady_state(ss_V, x_found[:n_params])
    else:
        x = mm.steady_state(ss_V, x_found[n_params:2*n_params])
    for i in range(len(states)):
        m.get(states[i]).set_state_value(x[i])

    log = data.load_ap_protocol().npview()
    t, v = log['time'], log['voltage']

    s = myokit.Simulation(m, p)
    s.set_fixed_form_protocol(t, v)
    s.set_tolerance(1e-8, 1e-8)
    s.set_max_step_size(0.1)

    # Update model parameters
    for i in range(n_params):
        if mutant in {1, 2}:
            s.set_constant('ikr.p'+str(i+1), x_found[i])
        else:
            s.set_constant('ikr.p'+str(i+1), x_found[i+n_params])

    d = s.run(p.characteristic_time(), log_interval=dt, log=d)

    signals2 = [d.time(), d['ikr.IKr'], d['membrane.V']]
    d = myokit.DataLog()
    d.set_time_key('time')
    d['time'] = signals2[0]
    d['current'] = signals2[1]
    d['voltage'] = signals2[2]

    # Filtered simulated data
    d = d.npview()
    ds.append(d)

d_WT = ds[0]
d_R56Q = ds[1]
d_RWT = ds[2]
d_RP56Q = ds[3]

a0.set_xlim([0, 8000])
a0.set_ylabel('Voltage (mV)', fontsize=9)
a0.plot(d_WT.time(), d_WT['voltage'], color='black')
a0.grid(True)
[label.set_visible(False) for label in a0.get_xticklabels()]
a1.axis([0, 8000, -14, 6])
a1.set_xlabel('Time (ms)', fontsize=9)
a1.set_ylabel('Current (nA)', fontsize=9)
a1.plot(d_WT.time(), d_WT['current'], c=(0.11765,0.23529,1.0))
a1.plot(d_R56Q.time(), d_R56Q['current'], c='red')
a1.legend(['WT', 'R56Q'], loc='upper right', fontsize=8, ncol=3)
a1.grid(True)
axins = zoomed_inset_axes(a1, 2.7, loc='lower left')
axins.plot(d_WT.time(), d_WT['current'], c=(0.11765,0.23529,1.0))
axins.plot(d_R56Q.time(), d_R56Q['current'], c='red')
x1, x2, y1, y2 = 2550, 3350, -0.2, 4.6 # specify the limits
axins.set_xlim(x1, x2) # apply the x-limits
axins.set_ylim(y1, y2) # apply the y-limits
axins.grid(True)
pl.yticks(visible=False)
pl.xticks(visible=False)
mark_inset(a1, axins, loc1=2, loc2=1, fc="none", ec="0.5")

axins = zoomed_inset_axes(a1, 3.2, loc='lower center')
axins.plot(d_WT.time(), d_WT['current'], c=(0.11765,0.23529,1.0))
axins.plot(d_R56Q.time(), d_R56Q['current'], c='red')
x1, x2, y1, y2 = 5020, 5920, -0.2, 3.6 # specify the limits
axins.set_xlim(x1, x2) # apply the x-limits
axins.set_ylim(y1, y2) # apply the y-limits
axins.grid(True)
pl.yticks(visible=False)
pl.xticks(visible=False)
mark_inset(a1, axins, loc1=2, loc2=1, fc="none", ec="0.5")

axins = zoomed_inset_axes(a1, 2.8, loc='lower right')
axins.plot(d_WT.time(), d_WT['current'], c=(0.11765,0.23529,1.0))
axins.plot(d_R56Q.time(), d_R56Q['current'], c='red')
x1, x2, y1, y2 = 6550, 7300, -0.2, 3.8 # specify the limits
axins.set_xlim(x1, x2) # apply the x-limits
axins.set_ylim(y1, y2) # apply the y-limits
axins.grid(True)
pl.yticks(visible=False)
pl.xticks(visible=False)
mark_inset(a1, axins, loc1=2, loc2=1, fc="none", ec="0.5")

if args.show == True:
    pl.show()
else:
    filename = 'AP/AP-R56Q-compare-model-' + model_str
    pl.savefig('PNG_figures/' + filename + '.png')
