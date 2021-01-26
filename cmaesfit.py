#!/usr/bin/env python3
#
from __future__ import division, print_function
import os
import sys
import pints
import numpy as np
import myokit
import argparse

# Check input arguments
parser = argparse.ArgumentParser(
    description='Fit a hERG model to electrophysiology data')
parser.add_argument('--model', type=int, default=2, metavar='N',
                    help='model number : 1 for CCOI, 2 for M10')
parser.add_argument('--mutant', type=int, default=1, metavar='N',
                    help='mutant number : 1 for WT, 2 for R56Q')
parser.add_argument('-r', '--repeats', type=int, default=20, metavar='N',
                    help='number of CMA-ES runs from different initial guesses')
parser.add_argument('-p', '--protocol', type=int, default=2, metavar='N',
                    help='which protocol is used to fit the data: 1 for sine wave, 2 for staircase #1, etc.')
args = parser.parse_args()

# Load project modules
sys.path.append(os.path.abspath(os.path.join('python')))
import priors
import cells
import transformation
import data
import model

# Fix seed
np.random.seed(1)

# Get string and protocols for each mutant
# Set Ek based on cell-specific estimate
if args.mutant == 1:
    mutant_str = 'WT'
else:
    mutant_str = 'R56Q'

protocol = args.protocol

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
elif args.protocol == 7:
    protocols = [2, 7]
    protocol_str = 'staircase1-conductance'
else:
    protocols = [7]
    protocol_str = 'conductance'

if args.mutant == 1:
    mutants = ['WT', 'WT-RPR']
    cs = [1, 2, 4]
else:
    mutants = ['R56Q', 'R56Q-RPR']
    cs = [1, 5, 6]

no_cells = len(cs)

# Get model string and params
if args.model == 1:
    model_str = 'CCOI'
    x_found = np.loadtxt('cmaesfits/parameter-sets/CCOI-params-' + mutant_str + '.txt', unpack=True)
elif args.model == 2:
    model_str = 'M10'
    x_found = np.loadtxt('cmaesfits/parameter-sets/M10-params-' + mutant_str + '.txt', unpack=True)
else:
    print('Invalid model')
    sys.exit()

trans = transformation.Transformation()
x_found = x_found
for i in range(no_cells):
    x_found = np.append(x_found, 0.07)
x_initial = trans.transform(no_cells=no_cells, parameters=x_found, which_model=args.model)
ek = cells.ek_computed()

filename = 'cmaesfits/' + mutant_str + '-model-' + model_str + '-fit-' + protocol_str + '-artefact'
print('Selected model ' + model_str)
print('Selected mutant ' + mutant_str)
print('Selected fitting protocol ' + protocol_str)
print('Storing results to ' + filename + '.txt')

# Define problems
models = []
for cell in cs:

    for mutant_str in mutants:

        # Create protocol
        if protocol in {4, 5}:
            p = data.load_myokit_protocol(protocol)
        else:
            p = data.load_protocol_values(protocol, mutant_str, cell)

        # Create forward models
        models.append(model.Model(
            p,
            cell=cell,
            cells=cs,
            EK=ek,
            mutant=mutant_str,
            which_model=args.model,
            analytical=(protocol in {4, 5}),
            conductance=(protocol == 7),
            sine_wave=(protocol == 1),
            staircase=(protocol in {2, 3}),
            step=(cell in {4, 5} and mutant_str in {'WT', 'WT-RPR'}),
            transformation=trans
        ))

# Load data, create single output problems
problems = []
sigma_noise = []
k = 0
for cell in cs:
    for mutant_str in mutants:

        log = data.load(cell, protocol, mutant_str)
        time = log.time()
        current = log['current'] / 1000 #change units from pA to nA
        
        debug = False
        if debug:
            test = models[k].simulate(x_initial, time)
            import matplotlib.pyplot as plt
            plt.figure()
            plt.plot(time, current)
            plt.plot(time, test)
            plt.grid(True)
            plt.show()

        # Estimate noise from the first 100 ms of data
        sigma_noise.append(np.std(current[:1000], ddof=1))
        problems.append(pints.SingleOutputProblem(models[k], time, current))
        k += 1
        del(log)

# Define log-likelihood

# Create log-posterior
log_likelihoods = []
for k, problem in enumerate(problems):
    print('Sigma noise', sigma_noise[k])
    log_likelihoods.append(pints.GaussianKnownSigmaLogLikelihood(problem, sigma_noise[k]))
log_likelihood = pints.SumOfIndependentLogPDFs(log_likelihoods)
log_prior = priors.LogPrior(no_cells=no_cells, which_model=args.model, transformation=trans)
log_posterior = pints.LogPosterior(log_likelihood, log_prior)
f = log_posterior

print('Score at initial parameters: ',
    f(x_initial))

#
# Run
#
b = myokit.Benchmarker()
repeats = args.repeats
params, scores = [], []
times = []
for i in range(repeats):
    print('Repeat ' + str(1 + i))

    # Choose random starting point
    if i < 18:
        fac = np.random.normal(1.0, 0.1, len(x_initial))
        q0 = fac*x_initial
    else:
        q0 = log_prior.sample()    # Search space

    # Create optimiser
    opt = pints.OptimisationController(
        f, q0, method=pints.CMAES)
    opt.set_log_to_file(filename + '-log-' + str(i) + '.txt')
    opt.set_max_iterations(None)
    opt.set_parallel(True)

    # Run optimisation
    try:
        with np.errstate(all='ignore'): # Tell numpy not to issue warnings
            b.reset()
            q, s = opt.run()            # Search space
            times.append(b.time())
            p = trans.detransform(q,args.model)    # Model space
            params.append(p)
            scores.append(-s)
    except ValueError:
        import traceback
        traceback.print_exc()

# Order from best to worst
order = np.argsort(scores)
scores = np.asarray(scores)[order]
params = np.asarray(params)[order]
times = np.asarray(times)[order]

# Show results
print('Best scores:')
for score in scores[:10]:
    print(-score)
print('Mean & std of score:')
print(-np.mean(scores))
print(np.std(scores))
print('Worst score:')
print(scores[-1])

# Extract best
obtained_score = scores[0]
obtained_parameters = params[0]

# Store results
print('Storing best result...')
with open(filename + '.txt', 'w') as f:
    for x in obtained_parameters:
        f.write(pints.strfloat(x) + '\n')

print('Storing all errors')
with open(filename + '-errors.txt', 'w') as f:
    for score in scores:
        f.write(pints.strfloat(-score) + '\n')

print('Storing all parameters')
for i, param in enumerate(params):
    with open(filename + '-parameters-' + str(1 + i) + '.txt', 'w') as f:
        for x in param:
            f.write(pints.strfloat(x) + '\n')

print('Storing all simulation times')
with open(filename + '-times.txt', 'w') as f:
    for time in times:
        f.write(pints.strfloat(time) + '\n')

