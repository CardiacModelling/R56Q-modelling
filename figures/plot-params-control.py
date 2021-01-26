import numpy as np
import argparse
import myokit
import matplotlib.pyplot as p
from matplotlib import *

p.rcdefaults()

p.rc('lines', markeredgewidth=2)
p.rc('lines', markersize=7)
p.rc('xtick.major', size=5)
p.rc('ytick.major', size=5) #changes size of y axis ticks on both sides
p.rc('xtick', direction='out')
p.rc('ytick', direction='out')

# Check input arguments
parser = argparse.ArgumentParser(
    description='Plot model parameters')
parser.add_argument('--model', type=int, default=2, metavar='N',
                    help='model number : 1 for CCOI, 2 for M10')
parser.add_argument('-p', '--protocol', type=int, default=2, metavar='N',
                    help='which protocol is used to fit the data: 1 for sine wave, 2 for staircase #1')
parser.add_argument("--show", action='store_true',
                    help="whether to show figures instead of saving them",
                    default=False)
parser.add_argument('--dpi', type=int, default=100, metavar='N',
                    help='what DPI to use for figures (suggested 300 for publication quality)')
args = parser.parse_args()

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

if args.model == 1:
    m = myokit.load_model('../model-and-protocols/CCOI-ikr-markov-voffset.mmt')
    model_str = 'CCOI'
elif args.model == 2:
    m = myokit.load_model('../model-and-protocols/M10-ikr-markov-voffset.mmt')
    model_str = 'M10'
else:
    print('Invalid model')
    sys.exit()   

n_kparams = int(m.get('misc.n_params').value())

colors = [(0.11765,0.23529,1.0), 'red']

filename = '../cmaesfits/WT-model-' + model_str + '-fit-' + protocol_str + '-artefact.txt'
params = np.loadtxt(filename, unpack=True)
filename = '../cmaesfits/R56Q-model-' + model_str + '-fit-' + protocol_str + '-artefact.txt'
params56Q = np.loadtxt(filename, unpack=True)

fig = p.figure(figsize=(7, 1.5), dpi=args.dpi, constrained_layout=True) #(8,6) seems to be the default
fig.set_facecolor('white') #this changes the background colour from the default grey/blue to white

grid = p.GridSpec(1, 1, figure=fig)

ax1 = p.subplot(grid[0, 0])
ax1.grid(True)
ax1.set_axisbelow(True)
ax1.semilogy(params[:n_kparams], marker='o', color=colors[0], linewidth=0, alpha=1, label='WT')
ax1.semilogy(params56Q[:n_kparams], marker='o', color=colors[1], linewidth=0, alpha=1, label='R56Q')
l = np.arange(0, n_kparams)
x = []
for i in range(n_kparams):
    x.append('p' + str(i+1))
ax1.set_xticks(l)
ax1.set_xticklabels(x)
ax1.legend(loc='lower right', ncol=2)

if args.show == True:
    p.show()
else:
    filename = 'Params-Control-model-' + model_str
    p.savefig('PNG_figures/' + filename + '.png')
