#!/usr/bin/env python3
#
# Pints Boundaries that limit the transition rates in the Beattie et al model.
#
from __future__ import division, print_function
import numpy as np
import pints


class LogPrior(pints.LogPrior):
    """
    Boundary constraints on the parameters
    """
    def __init__(self, no_cells=1, which_model=1, transformation=None):
        super(LogPrior, self).__init__()

        self.which_model = which_model

        # Conductance limits
        self.lower_conductance = 2e-2
        self.upper_conductance = 2e-1

        # Limits on p1-p8
        self.lower_alpha = 1e-7              # Kylie: 1e-7
        self.upper_alpha = 1e3               # Kylie: 1e3
        self.lower_beta  = 1e-7              # Kylie: 1e-7
        self.upper_beta  = 0.4               # Kylie: 0.4
        self.lower_rate  = 1.67e-5
        self.upper_rate  = 1000

        # Lower and upper bounds for all parameters
        self.lower = np.array([
            self.lower_alpha,
            self.lower_beta,
            self.lower_alpha,
            self.lower_beta,
            self.lower_alpha,
            self.lower_beta,
            self.lower_alpha,
            self.lower_beta,
            self.lower_alpha,
            self.lower_beta,
            self.lower_alpha,
            self.lower_beta
        ])
        self.upper = np.array([
            self.upper_alpha,
            self.upper_beta,
            self.upper_alpha,
            self.upper_beta,
            self.upper_alpha,
            self.upper_beta,
            self.upper_alpha,
            self.upper_beta,
            self.upper_alpha,
            self.upper_beta,
            self.upper_alpha,
            self.upper_beta
        ])

        self.minf = -float('inf')

        # Limits on maximum reaction rates
        self.rmin = 1.67e-5
        self.rmax = 1000

        # Voltages used to calculate maximum rates
        self.vmin = -120
        self.vmax =  60

        # Optional transformation
        self.transformation = transformation

        # Number of parameters
        n_kparams = 12

        self.no_cells = no_cells
        self.n_kparams = n_kparams
        self.n_params = 2*self.n_kparams + self.no_cells

        self.lower = np.append(self.lower, self.lower)
        self.upper = np.append(self.upper, self.upper)

        self.no_cells = no_cells
        for i in range(self.no_cells):
            self.lower = np.append(self.lower, self.lower_conductance)
            self.upper = np.append(self.upper, self.upper_conductance)

    def n_parameters(self):
        return self.n_params

    def __call__(self, parameters):

        debug = False

        # Transform parameters back to model space
        if self.transformation is not None:
            parameters = self.transformation.detransform(parameters, self.which_model)

        # Check parameter boundaries
        if np.any(parameters < self.lower):
            if debug: print('Lower')
            return self.minf
        if np.any(parameters > self.upper):
            if debug: print('Upper')
            return self.minf

        # Check maximum rate constants
        p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, \
        p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24 = parameters[:2*self.n_kparams]

        # Check positive signed rates
        r = p1 * np.exp(p2 * self.vmax)
        if r < self.rmin or r > self.rmax:
            if debug: print('r1')
            return self.minf
        r = p5 * np.exp(p6 * self.vmax)
        if r < self.rmin or r > self.rmax:
            if debug: print('r2')
            return self.minf
        r = p9 * np.exp(p10 * self.vmax)
        if r < self.rmin or r > self.rmax:
            if debug: print('r3')
            return self.minf

        # Check negative signed rates
        r = p3 * np.exp(-p4 * self.vmin)
        if r < self.rmin or r > self.rmax:
            if debug: print('r4')
            return self.minf
        r = p7 * np.exp(-p8 * self.vmin)
        if r < self.rmin or r > self.rmax:
            if debug: print('r5')
            return self.minf  
        r = p11 * np.exp(-p12 * self.vmin)
        if r < self.rmin or r > self.rmax:
            if debug: print('r6')
            return self.minf 

        # Check positive signed rates
        r = p13 * np.exp(p14 * self.vmax)
        if r < self.rmin or r > self.rmax:
            if debug: print('r7')
            return self.minf
        r = p17 * np.exp(p18 * self.vmax)
        if r < self.rmin or r > self.rmax:
            if debug: print('r8')
            return self.minf
        r = p21 * np.exp(p22 * self.vmax)
        if r < self.rmin or r > self.rmax:
            if debug: print('r9')
            return self.minf

        # Check negative signed rates
        r = p15 * np.exp(-p16 * self.vmin)
        if r < self.rmin or r > self.rmax:
            if debug: print('r10')
            return self.minf
        r = p19 * np.exp(-p20 * self.vmin)
        if r < self.rmin or r > self.rmax:
            if debug: print('r11')
            return self.minf  
        r = p23 * np.exp(-p24 * self.vmin)
        if r < self.rmin or r > self.rmax:
            if debug: print('r12')
            return self.minf    

        return True

    def _sample_partial(self, v):
        """
        Sample a pair of parameters - uniformly in the transformed space - that
        satisfy the maximum transition rate constraints.
        """
        for i in range(100):
            a = np.exp(np.random.uniform(
                np.log(self.lower_alpha), np.log(self.upper_alpha)))
            b = np.random.uniform(self.lower_beta, self.upper_beta)
            r = a * np.exp(b * v)
            if r >= self.rmin and r <= self.rmax:
                return a, b
        raise ValueError('Too many iterations')

    def sample(self, n=1):

        if n > 1:
            raise NotImplementedError

        p = np.zeros(self.n_params)

        # Sample forward rates
        p[0:2] = self._sample_partial(self.vmax)
        p[4:6] = self._sample_partial(self.vmax)
        p[8:10] = self._sample_partial(self.vmax)

        p[12:14] = self._sample_partial(self.vmax)
        p[16:18] = self._sample_partial(self.vmax)
        p[20:22] = self._sample_partial(self.vmax)

        # Sample backward rates
        p[2:4] = self._sample_partial(-self.vmin)
        p[6:8] = self._sample_partial(-self.vmin)
        p[10:12] = self._sample_partial(-self.vmin)

        p[14:16] = self._sample_partial(-self.vmin)
        p[18:20] = self._sample_partial(-self.vmin)
        p[22:24] = self._sample_partial(-self.vmin)

        for i in range(self.no_cells):
            p[2*self.n_kparams+i] = np.random.uniform(
                self.lower_conductance, self.upper_conductance)

        # Transform from model to search space, if required
        if self.transformation is not None:
            p = self.transformation.transform(no_cells=self.no_cells, parameters=p, which_model=self.which_model)

        # The Boundaries interface requires a matrix ``(n, n_parameters)``
        p.reshape(1, self.n_params)
        
        return p
        
