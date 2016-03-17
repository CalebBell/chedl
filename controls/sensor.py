# -*- coding: utf-8 -*-
'''Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, Caleb Bell <Caleb.Andrew.Bell@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''

from __future__ import division
import ctypes
import scipy.stats as stats






__random_distributions = ['anglit', 'arcsine', 'cauchy', 'cosine', 'expon',
                          'gilbrat', 'halfcauchy', 'halflogistic', 'halfnorm',
                          'hypsecant', 'kstwobign', 'laplace', 'levy',
                          'logistic', 'maxwell', 'norm', 'rayleigh',
                          'semicircular', 'uniform', 'wald']


class Sensor(object):
    '''
    Sensor package includes transmitter, variable reader, and settings.
    Can `output` the variable directly, or 4-20 mA or a-b V signal.
    `sensed variable` is the desired information, i.e. T for a thermocouple.
    `sensed signal

    The `output` max and min are different than the `variable` min and max;
    functional relationships normally allow the `variable` to be calculated
    in a much larger range than the `output` can be set. If the `output`
    range is not set, its range is set to those of the variable.
    '''
    t = 0 # Initialize at time = 0
    linear = False # Assume non-linear signal-variable relation

    def __init__(self, output_max=None, output_min=None, output_type='variable',
                 output_units=None, output_action='direct',
                 output_rounding_algorithm=None,
                 variable_min=None, variable_max=None, variable_units=None,
                 signal_min=None, signal_max=None, signal_units=None,
                 signal_rounding_algorithm=None, signal_error_distribution='norm',
                 deadtime=None, inaccuracy=None, offset=None):


         self.variable_min = self.variable_zero = variable_min
         self.variable_max = variable_max
         self.variable_span = self.variable_max - self.variable_min
         self.variable_units = variable_units

         self.signal_min = self.signal_zero = signal_min
         self.signal_max = signal_max
         self.signal_span = self.signal_max - self.signal_min
         self.signal_units = signal_units
         self.signal_rounding_algorithm = signal_rounding_algorithm
         self.signal_error_distribution = signal_error_distribution

         self.output_max = output_max
         self.output_min = output_min
         self.output_type = output_type
         self.output_units = output_units
         self.output_action = output_action
         self.output_rounding_algorithm = output_rounding_algorithm
         self.__initialize_outputs__()
         # Must add output units here too, based on output type

         self.deadtime = deadtime
         self.inaccuracy = inaccuracy
         self.offset = offset

    def __initialize_outputs__(self):
        if self.output_type == 'variable':
            self.output_min = self.output_zero = self.variable_min
            self.output_max = self.variable_max
            self.output_units = self.variable_units
        elif self.output_type == '4-20 mA':
            self.output_min = self.output_zero = 4.
            self.output_max = 20.
            self.output_span = self.output_max - self.output_min
            self.output_units = 'mA'

        elif self.output_type == '1-5 V':
            self.output_min = self.output_zero = 1.
            self.output_max = 5.
            self.output_span = self.output_max - self.output_min
            self.output_units = 'V'
        elif self.output_type == '0-5 V':
            self.output_min = self.output_zero = 0.
            self.output_max = 5.
            self.output_span = self.output_max - self.output_min
            self.output_units = 'V'
        else:
            raise ValueError(self.output_type)

    def output(self, signal):
        self.variable = self.signal_to_variable(signal)
        self.output = self.variable_to_output(self.variable)
        return self.output

        # Dummy function, to be overwridden by those in subclasses.
        return signal

    def variable_to_output(self, variable):
        if self.output_type == 'variable':
            output = variable
        elif self.output_action == 'direct':
            output = ((variable-self.variable_min)
            /(self.variable_max - self.variable_min)
            *(self.output_max - self.output_min)) + self.output_min
        elif self.output_action == 'reverse':
            output = -((variable-self.variable_min)
            /(self.variable_max - self.variable_min)
            *(self.output_max - self.output_min)) + self.output_max
        else:
            raise Exception(ValueError)
        if self.output_rounding_algorithm:
            output = self.rounding(output, output_rounding_algorithm, self.output_max, self.output_min)

        return output

    def rounding(self, unrounded, rounding_algorithm=None, rounded_max=None, rounded_min=None):
        if rounding_algorithm == 'int':
            output = int(unrounded)
        elif rounding_algorithm == 'int8':
            output = ctypes.c_int8(unrounded).value
        elif rounding_algorithm == 'int16':
            output = ctypes.c_int16(unrounded).value
        elif rounding_algorithm == 'int32':
            output = ctypes.c_int32(unrounded).value
        elif rounding_algorithm == 'float':
            output = ctypes.c_float(unrounded).value
        elif 'bit' in rounding_algorithm:
            bits = int(rounding_algorithm.replace('bit', ''))
            output = unrounded - unrounded % (rounded_max - rounded_min)/2.**bits
        else:
            raise Exception(ValueError)
        return output

    def signal_error(self, signal):
        if self.inaccuracy and self.signal_error_distribution:
            signal += signal*getattr(stats, dist).rvs(loc=inaccuracy)
        return signal
