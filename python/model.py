#!/usr/bin/env python3
#
# Pints ForwardModel that runs simulations with the Markov herg models.
#
from __future__ import division, print_function
import myokit
import myokit.lib.markov as markov
import numpy as np
import pints

from data import model_path


class Model(pints.ForwardModel):
    """
    Pints ForwardModel that runs simulations with M10 or CCOI model.

    Arguments:

        ``protocol``
            A myokit.Protocol or a tuple (times, voltage)
        ``EK``
            The reversal potential
        ``sine_wave``
            Set to True if sine-wave protocol is being used.
        ``start_steady``
            Start at steady state for 0, -40, or -80mV.
        ``analytical``
            Use an analytical simulation.

    """

    def __init__(
        self, protocol, cell, cells, mutant, EK=-88.0, which_model=1, analytical=False, conductance=False, sine_wave=False,\
        staircase=False, start_steady=True, step=False, transformation=None):

        # Load model
        self.which_model = which_model
        if self.which_model == 1:
            model = myokit.load_model(model_path('CCOI-ikr-markov-voffset.mmt'))
            self.states = [
                'ikr.c1',
                'ikr.c2',
                'ikr.i',
                'ikr.o'
            ]
        else:
            model = myokit.load_model(model_path('M10-ikr-markov-voffset.mmt'))
            self.states = [
                'ikr.c1',
                'ikr.c2',
                'ikr.i',
                'ikr.ic1',
                'ikr.ic2',
                'ikr.o'
            ]
        self.model = model

        n_kparams = int(self.model.get('misc.n_params').value())
        no_cells = len(cells)
        parameters = np.zeros(2*n_kparams + no_cells)

        for i in range(n_kparams):
            parameters[i] = self.model.get('ikr.p' + str(i+1)).value()
            parameters[i+n_kparams] = parameters[i]            

        for i in range(no_cells):
            parameters[2*n_kparams+i] = 0.1

        self.sine_wave = sine_wave
        self.staircase = staircase
        self.step = step
        self.start_steady = start_steady

        self.parameters = parameters
        self.n_params = len(self.parameters)
        LJP = self.model.get('misc.LJP').value()

        if step:
            self.ss_V = 0
        elif conductance:
            self.ss_V = 40
        else:
            self.ss_V = -80
        self.ss_V = self.ss_V - LJP

        # Set reversal potential
        print('EK:', EK)
        self.model.get('nernst.EK').set_rhs(EK)

        if mutant in {'WT', 'WT-RPR'}:
            Voffs = [1.7, 1.255, -1.47, 1.63, 0.0]
            all_cells = [1, 2, 3, 4, 5]
        else:
            Voffs = [0.37, 0.0, 0.0, 0.61, 0.36, -1.39]
            all_cells = [1, 2, 3, 4, 5, 6]

        cell_idx = cells.index(cell)
        all_cell_idx = all_cells.index(cell)
        self.model.get('ikr.p' + str(n_kparams+1)).set_rhs(parameters[2*n_kparams+cell_idx])
        self.model.get('voltage_clamp.V_off').set_rhs(Voffs[all_cell_idx])
        self.cell_idx = cell_idx

        # Start at steady-state for -80 mV for Claydon 37C data
        if self.start_steady:
            print('Updating model to steady-state for ' + str(self.ss_V) + ' mV')
            self.model.get('membrane.V').set_label('membrane_potential')
            mm = markov.LinearModel.from_component(self.model.get('ikr'))
            self.mm = mm
            # Update states
            if mutant in {'WT', 'R56Q'}:
                x = mm.steady_state(self.ss_V, self.parameters[:(n_kparams+1)])
            else:
                x = mm.steady_state(self.ss_V, self.parameters[n_kparams:(2*n_kparams+1)])
            for i in range(len(self.states)):
                self.model.get(self.states[i]).set_state_value(x[i])

        if self.sine_wave:
            self.model.get('membrane.V').set_rhs(
                'piecewise(engine.time >= 300.0 and engine.time < 900.0,'
                + ' -140 - misc.LJP'
                + ' + 0.1 * (engine.time - 300.0),'
                + 'engine.time >= 3599.9 and engine.time < 7100.1,'
                + ' - 30 - misc.LJP'
                + ' + 54 * sin(0.007 * (engine.time - 3100.1))'
                + ' + 26 * sin(0.037 * (engine.time - 3100.1))'
                + ' + 10 * sin(0.190 * (engine.time - 3100.1)),'
                + 'engine.time >= 7100.1 and engine.time < 7200.1,'
                + ' -70 - misc.LJP'
                + ' - 0.3 * (engine.time - 7100.0)'
                + ', engine.pace - misc.LJP)')

        if self.staircase:
            if self.step:
                self.model.get('membrane.V').set_rhs(
                    'piecewise(engine.time <= 1296.8,'
                    + ' 0 - misc.LJP,'
                    + 'engine.time >= 14410.1 and engine.time < 14510.0,'
                    + ' -70 - misc.LJP'
                    + ' - 0.4 * (engine.time - 14410.1)'
                    + ', engine.pace - misc.LJP)')
            else:
                self.model.get('membrane.V').set_rhs(
                   'piecewise(engine.time <= 1236.2,'
                   + ' -80 - misc.LJP,'
                   + 'engine.time >= 14410.1 and engine.time < 14510.0,'
                   + ' -70 - misc.LJP'
                   + ' - 0.4 * (engine.time - 14410.1)'
                   + ', engine.pace - misc.LJP)')

        # Create simulation
        self._analytical = analytical
        if not self._analytical:
            self.simulation = myokit.Simulation(self.model)
            # Add protocol
            if isinstance(protocol, myokit.Protocol):
                self.simulation.set_protocol(protocol)
            else:
                # Apply data-clamp
                times, voltage = protocol
                self.simulation.set_fixed_form_protocol(times, voltage)

                # Set max step size
                self.simulation.set_max_step_size(0.1)

            # Set solver tolerances
            self.simulation.set_tolerance(1e-8, 1e-8)
        else:
            if self.sine_wave:
                raise ValueError(
                    'Analytical simulation cannot be used with sine wave or staircase protocols.')
            elif not isinstance(protocol, myokit.Protocol):
                raise ValueError(
                    'Analytical simulation cannot be used with data clamp.')
            if not self.start_steady:
                self.model.get('membrane.V').set_label('membrane_potential')
                mm = markov.LinearModel.from_component(model.get('ikr'))
            self.simulation = markov.AnalyticalSimulation(mm, protocol)

        self.transformation = transformation
        # Set a maximum duration for each simulation.
        self._timeout = myokit.Timeout(60)
        self.n_kparams = n_kparams
        self.cell = cell
        self.cells = cells
        self.no_cells = no_cells
        self.mutant = mutant
        self.protocol = protocol

    def n_parameters(self):
        return self.n_params

    def set_tolerances(self, tol):
        self.simulation.set_tolerance(tol, tol)

    def simulate(self, parameters, times):

        if self.transformation is not None:
            parameters = self.transformation.detransform(parameters, self.which_model)

        if self._analytical:
            self.simulation = markov.AnalyticalSimulation(self.mm, self.protocol)

        # Update model parameters
        for i in range(self.n_kparams):
            if self.mutant in {'WT', 'R56Q'}:
                self.simulation.set_constant('ikr.p'+str(i+1), parameters[i])
            else:
                self.simulation.set_constant('ikr.p'+str(i+1), parameters[i+self.n_kparams])

        self.simulation.set_constant('ikr.p'+str(self.n_kparams+1), parameters[2*self.n_kparams+self.cell_idx])

        # Run
        self.simulation.reset()
        try:
            if self._analytical:
                d = self.simulation.run(
                    times[-1] + times[1],
                    log_times=times,
                    ).npview()
            else:
                d = self.simulation.run(
                    times[-1] + times[1],
                    log_times=times,
                    log=['ikr.IKr', 'membrane.V'],
                    progress=self._timeout,
                    ).npview()
        except myokit.SimulationError:
            return times * float('inf')
        except myokit.SimulationCancelledError:
            return times * float('inf')

        # Return
        return d['ikr.IKr']

