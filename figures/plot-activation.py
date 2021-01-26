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
                    help='mutant number : 1 for WT, 2 for R56Q, 3 for WT+RPR, 4 for R56Q+RPR')
parser.add_argument('-p', '--protocol', type=int, default=2, metavar='N',
                    help='which protocol is used to fit the data: 1 for sine wave, 2 for staircase #1')
parser.add_argument("--show", action='store_true',
                    help="whether to show figures instead of saving them",
                    default=False)
parser.add_argument('--dpi', type=int, default=100, metavar='N',
                    help='what DPI to use for figures (suggested 300 for publication quality)')
args = parser.parse_args()

pr_steps = [
    (16000, 1000), #-50
    (46000, 1000), #-35
    (76000, 1000), #-20
    (106000, 1000), #-5
    (136000, 1000), #10
    (166000, 1000), #25
    (196000, 1000), #40
]
pr_voltages = np.array([-50,-35,-20,-5,10,25,40])

# Get model
p = myokit.load_protocol('../model-and-protocols/activation-ramp.mmt')

current = 'ikr.IKr'

if args.mutant == 1:
    mutant_str = 'WT'
else:
    mutant_str = 'R56Q'

ek = cells.ek_computed()

print('Reversal potential ' + str(ek) + ' mV')

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

# Number of steps in pacing protocol
n_steps = 7

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

if args.mutant in {1, 3}:
    cs = [1, 2, 4]
else:
    cs = [1, 5, 6]
no_cells = len(cs)

RPRS = ['no-RPR', 'RPR']

exp_ss_vs, exp_ss_ins, v1s, g1s = [], [], [], []
exp_iv_vs, exp_iv_ins, v2s, g2s = [], [], [], []

for c in cs:

    cell = c
    if args.mutant == 1:
        mutant_str = 'WT'
    else:
        mutant_str = 'R56Q'

    fig = pl.figure(figsize=(5, 4), dpi=args.dpi)

    ax1 = fig.add_subplot(3,2,1)
    pl.title('Control')
    ax1.set_ylabel( 'Voltage (mV)', fontsize=9 )
    [label.set_visible(False) for label in ax1.get_xticklabels()]
    ax1.grid(True)
    ax1.set_xlim([0, 3000])

    ax4 = fig.add_subplot(3,2,2)
    pl.title('RPR')
    [label.set_visible(False) for label in ax4.get_xticklabels()]
    [label.set_visible(False) for label in ax4.get_yticklabels()]
    ax4.grid(True)
    ax4.set_xlim([0, 3000])

    ax2 = fig.add_subplot(3,2,3)
    ax2.set_ylabel( 'Predicted\ncurrent (nA)', fontsize=9 )
    [label.set_visible(False) for label in ax2.get_xticklabels()]
    ax2.grid(True)
    ax2.set_xlim([0, 3000])

    ax5 = fig.add_subplot(3,2,4)
    [label.set_visible(False) for label in ax5.get_xticklabels()]
    ax5.grid(True)
    ax5.set_xlim([0, 3000])

    ax3 = fig.add_subplot(3,2,5)
    ax3.set_xlabel( 'Time (ms)', fontsize=9 )
    ax3.set_ylabel( 'Experimental\ncurrent (nA)', fontsize=9 )
    ax3.grid(True)
    ax3.set_xlim([0, 3000])

    ax6 = fig.add_subplot(3,2,6)
    ax6.set_xlabel( 'Time (ms)', fontsize=9 )
    ax6.grid(True)
    ax6.set_xlim([0, 3000])

    for r in RPRS:

        # Set steady state potential
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

        print('Updating model to steady-state for ' + str(ss_V) + ' mV')

        mm = markov.LinearModel(m, states, parameters, current)

        if r == 'no-RPR':
            x = mm.steady_state(ss_V, x_found[:n_params])
        else:
            mutant_str = mutant_str + '-RPR'
            x = mm.steady_state(ss_V, x_found[n_params:2*n_params])
        print(x)
        for i in range(len(states)):
            m.get(states[i]).set_state_value(x[i])

        m.get('membrane.V').set_rhs('engine.pace - misc.LJP')

        s = myokit.Simulation(m, p)
        s.set_tolerance(1e-8, 1e-8)
        s.set_max_step_size(0.1)

        # Update model parameters
        for i in range(n_params):
            if r == 'no-RPR':
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

        # Create colormap for plotting
        if args.mutant == 1:
            cmap = matplotlib.cm.get_cmap('viridis')
            colour = (0.11765,0.23529,1.0)
            exp_colour = 'limegreen'
        else:
            cmap = matplotlib.cm.get_cmap('plasma')
            colour = 'red'
            exp_colour = 'orange'
        norm = matplotlib.colors.Normalize(0, n_steps)

        # Filtered simulated data
        d = d.npview()
        df = d.fold(3000)

        e = myokit.DataLog.load_csv('../data/Activation/activation-' + mutant_str + '-cell-' + str(cell) + '.csv').npview()

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

        # Filtered experimental data
        e = e.npview()
        e = e.fold(3000)

        # Compute activation curve
        v1, g1 = biomarkers.steady_state_activation(pr_steps, pr_voltages, log=d, normalise=True)
        v2, g2 = biomarkers.iv_curve_activation(pr_steps, pr_voltages, log=d, normalise=True)

        ng1 = g1/g1[-1]
        ng2 = g2/np.max(g2)

        import scipy.interpolate as sp

        interp = sp.interp1d(v1, g1, kind='cubic')
        interp2 = sp.interp1d(v2, g2, kind='cubic')
        xnew = np.linspace(-50, 40, 50)

        a = np.loadtxt('Biomarkers/exp/act_ss_exp_' + mutant_str + '_cell' + str(cell) + '.txt', unpack=True)
        b = np.loadtxt('Biomarkers/exp/act_iv_exp_' + mutant_str + '_cell' + str(cell) + '.txt', unpack=True)

        exp_ss_v = a[0] - LJP
        exp_ss_i = a[1] / 1000  # get the right units
        exp_ss_in = a[2]

        exp_iv_v = b[0] - LJP
        exp_iv_i = b[1] / 1000  # get the right units
        exp_iv_in = b[2]

        exp_ss_vs.append(exp_ss_v)
        exp_ss_ins.append(exp_ss_in)
        v1s.append(v1)
        g1s.append(g1)

        exp_iv_vs.append(exp_iv_v)
        exp_iv_ins.append(exp_iv_in)
        v2s.append(v2)
        g2s.append(g2)

        if r == 'no-RPR':
            for k in range(n_steps):
                ax1.plot(df.time(), df['voltage',k], linewidth=1, color=cmap(norm(k)))
                ax2.plot(df.time(), df['current',k], linewidth=1, color=cmap(norm(k)))
                ax3.plot(e.time(), e['current',k], linewidth=1, color=cmap(norm(k)))
        else:
            for k in range(n_steps):
                ax4.plot(df.time(), df['voltage',k], linewidth=1, color=cmap(norm(k)))
                ax5.plot(df.time(), df['current',k], linewidth=1, color=cmap(norm(k)))
                ax6.plot(e.time(), e['current',k], linewidth=1, color=cmap(norm(k)))

        pl.tight_layout()

    if args.show == True:
        pl.show()
    else:
        filename = 'Activation/Activation-' + mutant_str + '-model-' + model_str + '-fit-' + protocol_str + '-cell-' + str(cell)
        pl.savefig('PNG_figures/' + filename + '.png')


defaults = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
if args.mutant == 2:
    defaults = defaults[3:]

fig = pl.figure(figsize=(5, 3), dpi=args.dpi)

ax1 = fig.add_subplot(2, 3, 1)
[label.set_visible(False) for label in ax1.get_xticklabels()]
ax1.set_ylabel('Normalised\ncurrent (control)', fontsize=9)
pl.scatter(exp_ss_vs[0], exp_ss_ins[0], s=50, facecolors='none', edgecolors=defaults[0], marker='o', linewidth=1)
pl.plot(v1s[0] - LJP, g1s[0], c=defaults[0], linestyle='--', linewidth=1)
ax1.grid(True)
ax1.text(0.8, 0.15, 'Cell 1',
        horizontalalignment='center',
        transform=ax1.transAxes)
ax2 = fig.add_subplot(2, 3, 2)
[label.set_visible(False) for label in ax2.get_xticklabels()]
[label.set_visible(False) for label in ax2.get_yticklabels()]
pl.scatter(exp_ss_vs[2], exp_ss_ins[2], s=50, facecolors='none', edgecolors=defaults[1], marker='o', linewidth=1)
pl.plot(v1s[2] - LJP, g1s[2], c=defaults[1], linestyle='--', linewidth=1)
ax2.grid(True)
ax2.text(0.8, 0.15, 'Cell 2',
        horizontalalignment='center',
        transform=ax2.transAxes)
ax3 = fig.add_subplot(2, 3, 3)
[label.set_visible(False) for label in ax3.get_xticklabels()]
[label.set_visible(False) for label in ax3.get_yticklabels()]
pl.scatter(exp_ss_vs[4], exp_ss_ins[4], s=50, facecolors='none', edgecolors=defaults[2], marker='o', linewidth=1)
pl.plot(v1s[4] - LJP, g1s[4], c=defaults[2], linestyle='--', linewidth=1)
ax3.grid(True)
ax3.text(0.8, 0.15, 'Cell 3',
        horizontalalignment='center',
        transform=ax3.transAxes)
ax4 = fig.add_subplot(2, 3, 4)
ax4.set_xlabel('Voltage (mV)', fontsize=9)
ax4.set_ylabel('Normalised\ncurrent (RPR)', fontsize=9)
pl.scatter(exp_ss_vs[1], exp_ss_ins[1], s=50, facecolors='none', edgecolors=defaults[0], marker='o', linewidth=1)
pl.plot(v1s[1] - LJP, g1s[1], c=defaults[0], linestyle='--', linewidth=1)
ax4.grid(True)
ax4.text(0.8, 0.15, 'Cell 1',
        horizontalalignment='center',
        transform=ax4.transAxes)
ax5 = fig.add_subplot(2, 3, 5)
[label.set_visible(False) for label in ax5.get_yticklabels()]
ax5.set_xlabel('Voltage (mV)', fontsize=9)
pl.scatter(exp_ss_vs[3], exp_ss_ins[3], s=50, facecolors='none', edgecolors=defaults[1], marker='o', linewidth=1)
pl.plot(v1s[3] - LJP, g1s[3], c=defaults[1], linestyle='--', linewidth=1)
ax5.grid(True)
ax5.text(0.8, 0.15, 'Cell 2',
        horizontalalignment='center',
        transform=ax5.transAxes)
ax6 = fig.add_subplot(2, 3, 6)
[label.set_visible(False) for label in ax6.get_yticklabels()]
ax6.set_xlabel('Voltage (mV)', fontsize=9)
pl.scatter(exp_ss_vs[5], exp_ss_ins[5], s=50, facecolors='none', edgecolors=defaults[2], marker='o', linewidth=1)
pl.plot(v1s[5] - LJP, g1s[5], c=defaults[2], linestyle='--', linewidth=1)
ax6.grid(True)
ax6.text(0.8, 0.15, 'Cell 3',
        horizontalalignment='center',
        transform=ax6.transAxes)

pl.tight_layout()

if args.show == True:
    pl.show()
else:
    filename = 'Biomarkers/Activation-SS-biomarker-' + mutant_str + '-model-' + model_str + '-fit-' + protocol_str
    pl.savefig('PNG_figures/' + filename + '.png')

fig = pl.figure(figsize=(5, 3), dpi=args.dpi)

ax1 = fig.add_subplot(2, 3, 1)
[label.set_visible(False) for label in ax1.get_xticklabels()]
ax1.set_ylabel('Normalised\ncurrent (control)', fontsize=9)
pl.scatter(exp_iv_vs[0], exp_iv_ins[0], s=50, facecolors='none', edgecolors=defaults[0], marker='o', linewidth=1)
pl.plot(v2s[0] - LJP, g2s[0], c=defaults[0], linestyle='--', linewidth=1)
ax1.grid(True)
ax1.text(0.8, 0.05, 'Cell 1',
        horizontalalignment='center',
        transform=ax1.transAxes)
ax2 = fig.add_subplot(2, 3, 2)
[label.set_visible(False) for label in ax2.get_xticklabels()]
[label.set_visible(False) for label in ax2.get_yticklabels()]
pl.scatter(exp_iv_vs[2], exp_iv_ins[2], s=50, facecolors='none', edgecolors=defaults[1], marker='o', linewidth=1)
pl.plot(v2s[2] - LJP, g2s[2], c=defaults[1], linestyle='--', linewidth=1)
ax2.grid(True)
ax2.text(0.8, 0.05, 'Cell 2',
        horizontalalignment='center',
        transform=ax2.transAxes)
ax3 = fig.add_subplot(2, 3, 3)
[label.set_visible(False) for label in ax3.get_xticklabels()]
[label.set_visible(False) for label in ax3.get_yticklabels()]
pl.scatter(exp_iv_vs[4], exp_iv_ins[4], s=50, facecolors='none', edgecolors=defaults[2], marker='o', linewidth=1)
pl.plot(v2s[4] - LJP, g2s[4], c=defaults[2], linestyle='--', linewidth=1)
ax3.grid(True)
ax3.text(0.8, 0.05, 'Cell 3',
        horizontalalignment='center',
        transform=ax3.transAxes)
ax4 = fig.add_subplot(2, 3, 4)
ax4.set_xlabel('Voltage (mV)', fontsize=9)
ax4.set_ylabel('Normalised\ncurrent (RPR)', fontsize=9)
pl.scatter(exp_iv_vs[1], exp_iv_ins[1], s=50, facecolors='none', edgecolors=defaults[0], marker='o', linewidth=1)
pl.plot(v2s[1] - LJP, g2s[1], c=defaults[0], linestyle='--', linewidth=1)
ax4.grid(True)
ax4.text(0.8, 0.05, 'Cell 1',
        horizontalalignment='center',
        transform=ax4.transAxes)
ax5 = fig.add_subplot(2, 3, 5)
[label.set_visible(False) for label in ax5.get_yticklabels()]
ax5.set_xlabel('Voltage (mV)', fontsize=9)
pl.scatter(exp_iv_vs[3], exp_iv_ins[3], s=50, facecolors='none', edgecolors=defaults[1], marker='o', linewidth=1)
pl.plot(v2s[3] - LJP, g2s[3], c=defaults[1], linestyle='--', linewidth=1)
ax5.grid(True)
ax5.text(0.8, 0.05, 'Cell 2',
        horizontalalignment='center',
        transform=ax5.transAxes)
ax6 = fig.add_subplot(2, 3, 6)
[label.set_visible(False) for label in ax6.get_yticklabels()]
ax6.set_xlabel('Voltage (mV)', fontsize=9)
pl.scatter(exp_iv_vs[5], exp_iv_ins[5], s=50, facecolors='none', edgecolors=defaults[2], marker='o', linewidth=1)
pl.plot(v2s[5] - LJP, g2s[5], c=defaults[2], linestyle='--', linewidth=1)
ax6.grid(True)
ax6.text(0.8, 0.05, 'Cell 3',
        horizontalalignment='center',
        transform=ax6.transAxes)

pl.tight_layout()

if args.show == True:
    pl.show()
else:
    filename = 'Biomarkers/Activation-IV-biomarker-' + mutant_str + '-model-' + model_str + '-fit-' + protocol_str
    pl.savefig('PNG_figures/' + filename + '.png')
