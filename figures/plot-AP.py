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

pl.rcParams.update({'figure.max_open_warning': 0})
pl.rcParams['axes.axisbelow'] = True

# Load project modules
sys.path.append(os.path.abspath(os.path.join('../', 'python')))
import cells
import data
import biomarkers

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
p = myokit.load_protocol('../model-and-protocols/RPR-protocol.mmt')

current = 'ikr.IKr'

protocol = args.protocol

if args.mutant == 1:
    mutant_str = 'WT'
else:
    mutant_str = 'R56Q'

ek = cells.ek_computed()

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
LJP = m.get('misc.LJP').value()
m.get('membrane.V').set_label('membrane_potential')

x_found = np.loadtxt('../cmaesfits/' + mutant_str + '-model-' + model_str + '-fit-' + protocol_str + '-artefact.txt', unpack=True)

if args.mutant == 1:
    cs = [1, 2, 4]
else:
    cs = [1, 5, 6]
no_cells = len(cs)

RPRS = ['no-RPR', 'RPR']

for c in cs:

    cell = c
    if args.mutant == 1:
        mutant_str = 'WT'
    else:
        mutant_str = 'R56Q'

    fig = pl.figure(c, figsize=(5, 4), dpi=args.dpi)
    ax1 = fig.add_subplot(3,1,1)
    ax1.set_ylabel('Voltage (mV)', fontsize=9)
    [label.set_visible(False) for label in ax1.get_xticklabels()]
    ax1.grid(True)
    ax1.set_xlim([0, 8000])
    ax2 = fig.add_subplot(3,1,2)
    ax2.set_ylabel('Current (nA)', fontsize=9)
    [label.set_visible(False) for label in ax2.get_xticklabels()]
    ax2.grid(True)
    ax2.set_xlim([0, 8000])
    ax3 = fig.add_subplot(3,1,3)
    ax3.set_xlabel('Time (ms)', fontsize=9)
    ax3.set_ylabel('Current (nA)', fontsize=9)
    ax3.grid(True)
    ax3.set_xlim([0, 8000])

    for r in RPRS:

        # Set steady state potential
        if cell in {4, 5} and mutant_str in {'WT', 'WT-RPR'}:
            ss_V = 0
        else:
            ss_V = -80

        ss_V = ss_V - LJP

        parameters = []
        for i in range(n_params):
            parameters.append('ikr.p'+str(i+1))

        d = [
            'engine.time',
            'membrane.V',
            'ikr.IKr'
            ]

        # Run simulation
        m.get('nernst.EK').set_rhs(ek)

        if mutant_str in {'WT', 'WT-RPR'}:
            Voffs = [1.7, 1.255, -1.47, 1.63, 0.0]
            all_cells = [1, 2, 3, 4, 5]
        else:
            Voffs = [0.37, 0.0, 0.0, 0.61, 0.36, -1.39]
            all_cells = [1, 2, 3, 4, 5, 6]

        cell_idx = cs.index(cell)
        all_cell_idx = all_cells.index(cell)
        g = x_found[2*n_params+cell_idx]
        Voff = Voffs[all_cell_idx]
        m.get('ikr.p' + str(n_params+1)).set_rhs(g)
        m.get('voltage_clamp.V_off').set_rhs(Voff)

        mm = markov.LinearModel(m, states, parameters, current)

        if r == 'no-RPR':
            x = mm.steady_state(ss_V, x_found[:n_params])
        else:
            mutant_str = mutant_str + '-RPR'
            x = mm.steady_state(ss_V, x_found[n_params:2*n_params])
        for i in range(len(states)):
            m.get(states[i]).set_state_value(x[i])

        log = data.load_RPR_protocol(cell, mutant_str).npview()
        t, v = log['time'], log['voltage']
        v -= LJP

        s = myokit.Simulation(m, p)
        s.set_fixed_form_protocol(t, v)
        s.set_tolerance(1e-8, 1e-8)
        s.set_max_step_size(0.1)

        # Update model parameters
        for i in range(n_params):
            if r == 'no-RPR':
                s.set_constant('ikr.p'+str(i+1), x_found[i])
            else:
                s.set_constant('ikr.p'+str(i+1), x_found[i+n_params])

        d = s.run(p.characteristic_time(), log_interval=dt, log=d)

        e = myokit.DataLog.load_csv('../data/ALL/all-' + mutant_str + '-cell-' + str(cell) + '.csv').npview()

        # Apply capacitance filtering for experiment and simulated data
        signals = [e.time(), e['current']]
        voltage = 'voltage' in e
        if voltage:
            signals.append(e['voltage'])
        signals = data.capacitance(p, dt, *signals)

        e = myokit.DataLog()
        e.set_time_key('time')
        e['time'] = signals[0]
        e['current'] = signals[1] / 1000
        if voltage:
            e['voltage'] = signals[2]

        signals2 = [d.time(), d['ikr.IKr'], d['membrane.V']]
        d = myokit.DataLog()
        d.set_time_key('time')
        d['time'] = signals2[0]
        d['current'] = signals2[1]
        d['voltage'] = signals2[2]

        # Create colormap for plotting
        if args.mutant == 1:
            colors = [(0.11765,0.23529,1.0), 'limegreen']
        else:
            colors = ['red', 'orange']

        # Filtered simulated data
        d = d.npview()

        # Filtered experimental data
        e = e.npview()

        d_staircase1 = d.trim_right(15400)
        e_staircase1 = e.trim_right(15400)

        d_activation = d.trim(15401,36401)
        e_activation = e.trim(15401,36401)

        d_factivation = d.trim(15401,36401,adjust=True)
        e_factivation = e.trim(15401,36401,adjust=True)

        d_ffactivation = d_factivation.fold(3000)
        e_ffactivation = e_factivation.fold(3000)

        d_inactivation = d.trim(36401,50402)
        e_inactivation = e.trim(36401,50401)

        d_finactivation = d.trim(36401,50402,adjust=True)
        e_finactivation = e.trim(36401,50401,adjust=True)

        d_ffinactivation = d_finactivation.fold(2000)
        e_ffinactivation = e_finactivation.fold(2000)

        d_ap = d.trim(50401,59225)
        e_ap = e.trim(50401,59225)

        d_sw = d.trim(59225,67225)
        e_sw = e.trim(59225,67225)

        d_validation = d.trim_left(15400)

        if r == 'no-RPR':
            RPR_str = 'control'
        else:
            RPR_str = 'RPR'

        if r == 'no-RPR':
            ax1.plot(d_ap.time() - d_ap.time()[0], d_ap['voltage'] + Voff, linewidth=1, color='black')
            ax2.plot(e_ap.time() - e_ap.time()[0], e_ap['current'], linewidth=1, color='silver', label='Expt. (' + RPR_str + ')')
            ax2.plot(d_ap.time() - d_ap.time()[0], d_ap['current'], linewidth=1, alpha=1, color=colors[0], label='Model')
            if mutant_str in {'WT', 'WT-RPR'}:
                ax2.legend(loc='lower center', fontsize=8, ncol=2)
            else:
                ax2.legend(loc='upper right', fontsize=8, ncol=2)
        else:
            ax3.plot(e_ap.time() - e_ap.time()[0], e_ap['current'], linewidth=1, color='silver', label='Expt. (' + RPR_str + ')')
            ax3.plot(d_ap.time() - d_ap.time()[0], d_ap['current'], linewidth=1, alpha=1, color=colors[1], label='Model')
            if mutant_str in {'WT', 'WT-RPR'}:
                ax3.legend(loc='lower center', fontsize=8, ncol=2) 
            else:
                ax3.legend(loc='upper right', fontsize=8, ncol=2)
        pl.tight_layout()

        if args.mutant == 1:
            mutant_str = 'WT'
        else:
            mutant_str = 'R56Q'        

    if args.show == True:
        pl.show()
    else:
        filename = 'AP-protocol-' + mutant_str + '-model-' + model_str + '-prediction-' + protocol_str + '-cell-' + str(cell)
        pl.savefig('PNG_figures/' + filename + '.png')
