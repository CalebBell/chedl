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
from math import log, exp, sqrt
from numpy import roots, imag, real
from controls.sensor import Sensor

# Coefs checked with original, from a NIST document
coefs_low = [-2.13534729, 3.18324720, -1.80143597, 0.71727204, 0.50344027, -0.61893395, -0.05332322, 0.28021362, 0.10715224, -0.29302865, 0.04459872, 0.11868632, -0.05248134]
coefs_high = [2.78157254, 1.64650916, -0.13714390, -0.00649767, -0.00234444, 0.00511868, 0.00187982, -0.00204472, -0.00046122, 0.00045724]
invcoefs_low = [0.183324722, 0.240975303, 0.209108771, 0.190439972, 0.142648498, 0.077993465, 0.012475611, -0.032267127, -0.075291522, -0.056470670, 0.076201285, 0.123893204, -0.029201193, -0.091173542, 0.001317696, 0.026025526]
invcoefs_high = [439.932854, 472.418020, 37.684494, 7.472018, 2.920828, 0.005184, -0.963864, -0.188732, 0.191203, 0.049025]

def RTD_resistance_ITS90(T):
    r'''Calculates the resisistance of a Resistance Temperature Detector,
    relative to its value at 273.15 K. Two ranges are used, one below 273.16 K
    until 13.8033 K and one above it, up to 961.78 degrees Celcius.

    .. math::
        \ln W_r(T_{90}) = A_0 + \sum_{i=1}^{12} A_i
        \left[\frac{\ln(T_{90}/273.16K) + 1.5}{1.5}\right]^i

        W_r(T_{90}) = C_0 + \sum_{i=1}^9 C_i\left[\frac{T_{90}/K
        - 754.15}{481}\right]^i

    Parameters
    ----------
    T : float
        Temperature, [K]

    Returns
    -------
    w : float
        Resistance of RTD relative to its value at 273.15 K, [-]

    Notes
    -----
    No warnings if outside the suggested range.

    Examples
    --------
    >>> RTD_resistance_ITS90(150.), RTD_resistance_ITS90(350.)
    (0.4984000571593147, 1.302905074678503)

    References
    ----------
    .. [1] Preston-Thomas, H. "The International Temperature Scale of 1990
       (ITS-90)." Metrologia 27, no. 1 (1990): 3.
       doi:10.1088/0026-1394/27/1/002.
    '''
    if T <= 273.16:
        tot = [coefs_low[i]*((log(T/273.16) + 1.5)/1.5)**i for i in range(1,len(coefs_low))]
        w = exp(coefs_low[0] + sum(tot))
    else:
        tot = [coefs_high[i]*((T-754.15)/481.)**i for i in range(1,len(coefs_high))]
        w = coefs_high[0] + sum(tot)
    return w


def RTD_temperature_ITS90(w):
    r'''Calculates the temperature of a Resistance Temperature Detector,
    given its resistance relative to its value at 273.15 K. Two ranges are
    used, one below 273.16 K until 13.8033 K and one above it, up to 961.78
    degrees Celcius.

    .. math::
        T_{90}/273.16K = B_0 + \sum_{i=1}^{15} B_i
        \left[\frac{W_r(T_{90})^{1/6} - 0.65}{0.35}\right]^i

        T_{90}/K - 273.15 = D_0 + \sum_{i=1}^9 D_i\left[\frac{W_r(T_{90})-
        2.64}{1.64}\right]^i

    Parameters
    ----------
    w : float
        Resistance of RTD relative to its value at 273.15 K, [-]

    Returns
    -------
    T : float
        Temperature, [K]

    Notes
    -----
    No warnings if outside the suggested range.
    The first funcntino is equivalent to the normal function to within 0.1 mK,
    and the second to 0.013 mK.

    Examples
    --------
    >>> RTD_temperature_ITS90(0.4984000571593147), RTD_temperature_ITS90(1.302905074678503)
    (150.00061624978588, 350.0000249661761)

    References
    ----------
    .. [1] Preston-Thomas, H. "The International Temperature Scale of 1990
       (ITS-90)." Metrologia 27, no. 1 (1990): 3.
       doi:10.1088/0026-1394/27/1/002.
    '''
    if w < 1.0000599918: # low
        tot = [invcoefs_low[i]*((w**(1/6.) -0.65)/0.35)**i for i in range(1,len(invcoefs_low))]
        T = (invcoefs_low[0] + sum(tot))*273.16
    else:
        tot = [invcoefs_high[i]*((w - 2.64)/1.64)**i for i in range(1,len(invcoefs_high))]
        T = (invcoefs_high[0] + sum(tot))+273.15
    return T



def Callendar_Van_Dusen_resistance(T, A=3.9083E-3, B=-5.775E-7, C=-4.183E-12):
    r'''Calculates the resisistance of a Resistance Temperature Detector,
    relative to its value at 273.15 K. Two ranges are used, one below 273.16 K
    until 73.15 K and one above it, up to 850 degrees Celcius.

    .. math::
        W_r = 1 + At + Bt^2 + Ct^3(t-100)

        W_r = 1 + At + Bt^2

        t = T - 273.15

    Parameters
    ----------
    T : float
        Temperature, [K]
    A, B, C : floats, optional
        Coefficients for the Callendar Van Dusen equation, [various]

    Returns
    -------
    w : float
        Resistance of RTD relative to its value at 273.15 K, [-]

    Notes
    -----
    No warnings if outside the suggested range.

    A = 3.9083E-3
    B = -5.775E-7,
    C = -4.183E-12

    Examples
    --------
    >>> Callendar_Van_Dusen_resistance(150.), Callendar_Van_Dusen_resistance(350.)
    (0.508191171034818, 1.2969421847562501)

    References
    ----------
    .. [1] IEC 60751 Standard, as given in Omega RTD Specifications, Technical
       Reference at http://www.omega.com/Temperature/pdf/RTDSpecs_Ref.pdf
    '''
    t = T - 273.15
    if t < 0:
        w = 1 + A*t + B*t**2 + C*t**3*(t-100)
    else:
        w = (1 + A*t + B*t**2)
    return w

#print RTD_resistance_ITS90(273.16), RTD_resistance_ITS90(200.)
#print [Callendar_Van_Dusen_temperature(150.), Callendar_Van_Dusen_temperature(350.)], 'hi'

def Callendar_Van_Dusen_temperature(w, A=3.9083E-3, B=-5.775E-7, C=-4.183E-12):
    r'''Calculates the temperature of a Resistance Temperature Detector,
    from its resistance relative to its value at 273.15 K. Two ranges are used,
    one below 273.16 K until 73.15 K and one above it, up to 850 degrees
    Celcius. The following equations are solved by a root-finding algorithm,
    or the quadratic equation is necessary.

    .. math::
        W_r = 1 + At + Bt^2 + Ct^3(t-100)

        W_r = 1 + At + Bt^2

        t = T - 273.15

    Parameters
    ----------
    w : float
        Resistance of RTD relative to its value at 273.15 K, [-]
    A, B, C : floats, optional
        Coefficients for the Callendar Van Dusen equation, [various]

    Returns
    -------
    T : float
        Temperature, [K]

    Notes
    -----
    No warnings if outside the suggested range.

    A = 3.9083E-3
    B = -5.775E-7,
    C = -4.183E-12

    Examples
    --------
    >>> Callendar_Van_Dusen_temperature(0.508191171034818), Callendar_Van_Dusen_temperature(1.2969421847562501)
    (150.00000000000006, 349.9999999999998)

    References
    ----------
    .. [1] IEC 60751 Standard, as given in Omega RTD Specifications, Technical
       Reference at http://www.omega.com/Temperature/pdf/RTDSpecs_Ref.pdf
    '''
    # Try the case of > 273.15 K
    # IEC 60751 standard says neglect A for positive temperatures
    t = (-A + (A**2 - 4*B*(1-w))**0.5)/(2*B)
    if t > 0:
        T = t +273.15
    else:
        ts = roots([C, -100*C, B, A, 1-w])
        T = [real(i)+273.15 for i in ts if not imag(i)][-1]
    return T
#


class RTD(Sensor):
    '''RTD Class'''
    def __init__(self,
                 variable_units='K', signal_units='ohm',
                 deadtime=0, offset=0, R0=100., method='IEC 60751'):
        self.R0 = R0
        self.method = method

        if self.method == 'ITS-90':
            variable_min = 13.81
            variable_max = 630.75 + 273.15
            signal_min = RTD_resistance_ITS90(variable_min)
            signal_max = RTD_resistance_ITS90(variable_max)
            self.signal_to_variable_function = RTD_temperature_ITS90
            self.variable_to_signal_function = RTD_resistance_ITS90
        elif 'IEC 60751' in self.method:
            variable_min = 73.15
            variable_max = 850 + 273.15
            signal_min = Callendar_Van_Dusen_resistance(variable_min)
            signal_max = Callendar_Van_Dusen_resistance(variable_max)
            self.signal_to_variable_function = Callendar_Van_Dusen_temperature
            self.variable_to_signal_function = Callendar_Van_Dusen_resistance
        else:
            raise Exception(ValueError)

        Sensor.__init__(self, variable_min=variable_min, variable_max=variable_max,
                 variable_units=variable_units, signal_min=signal_min*R0,
                 signal_max=signal_max*R0, signal_units=signal_units,
                 deadtime=deadtime, offset=offset)


    def signal_to_variable(self, signal):
        # Signal as sensed by the R0 resistor

        T = self.signal_to_variable_function(signal/self.R0)
        return T

