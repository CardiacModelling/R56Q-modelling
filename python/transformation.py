#!/usr/bin/env python3
#
from __future__ import division, print_function
import numpy as np


class Transformation(object):
    """
    Transforms from model to search space (and back).
    """

    def transform(self, parameters, which_model, no_cells):
        """
        Transform from model into search space.
        """
        x = np.array([
            np.log(parameters[0]),
            parameters[1],
            np.log(parameters[2]),
            parameters[3],
            np.log(parameters[4]),
            parameters[5],
            np.log(parameters[6]),
            parameters[7],
            np.log(parameters[8]),
            parameters[9],
            np.log(parameters[10]),
            parameters[11],
            np.log(parameters[12]),
            parameters[13],
            np.log(parameters[14]),
            parameters[15],
            np.log(parameters[16]),
            parameters[17],
            np.log(parameters[18]),
            parameters[19],
            np.log(parameters[20]),
            parameters[21],
            np.log(parameters[22]),
            parameters[23]
        ])   

        self.n_params = len(x)

        self.no_cells = no_cells
        for i in range(self.no_cells):
            x = np.append(x, parameters[self.n_params+i])

        return x

    def detransform(self, transformed_parameters, which_model, noise=False):
        """
        Transform back from search space to model space.
        """
        x = np.array([
            np.exp(transformed_parameters[0]),
            transformed_parameters[1],
            np.exp(transformed_parameters[2]),
            transformed_parameters[3],
            np.exp(transformed_parameters[4]),
            transformed_parameters[5],
            np.exp(transformed_parameters[6]),
            transformed_parameters[7],
            np.exp(transformed_parameters[8]),
            transformed_parameters[9],
            np.exp(transformed_parameters[10]),
            transformed_parameters[11],
            np.exp(transformed_parameters[12]),
            transformed_parameters[13],
            np.exp(transformed_parameters[14]),
            transformed_parameters[15],
            np.exp(transformed_parameters[16]),
            transformed_parameters[17],
            np.exp(transformed_parameters[18]),
            transformed_parameters[19],
            np.exp(transformed_parameters[20]),
            transformed_parameters[21],
            np.exp(transformed_parameters[22]),
            transformed_parameters[23]
        ])

        for i in range(self.no_cells):
            x = np.append(x, transformed_parameters[self.n_params+i])
        if noise:
            x = np.append(x, transformed_parameters[self.n_params+self.no_cells])

        return x

