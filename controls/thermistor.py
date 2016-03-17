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
# Thermistors are calibrated for specific resistances ONLY!

def Steinhart_Hart_temperature(R, A, B, C):
    r'''Calculates the temperature of a Thermistor, using its absolute
    resistance and Steinhard-Hart coefficients according to [1]_. Coefficients
    are calibrated for every single thermistor separately, due to large
    variations in production methods.

    .. math::
        \frac{1}{T} = A + B(\ln R) + C(\ln R)^3

    Parameters
    ----------
    R : float
        Resistance of thermistor at T [ohm]
    A, B, C : floats
        Coefficients for Steinhart-Hart function, [-]

    Returns
    -------
    T : float
        Temperature, [K]

    Notes
    -----
    No warnings produced if outside the suggested range.
    Thermistors are normally valid only for short ranges under 200 K.

    Examples
    --------
    Coefficients for Omega 44004 Thermistor, 2252 ohm at STP:

    >>> Steinhart_Hart_temperature(2252, A=1.468E-3, B=2.383E-4, C=1.007E-7)
    298.1604531883543

    References
    ----------
    .. [1] Steinhart, John S., and Stanley R. Hart. "Calibration Curves for
       Thermistors." Deep Sea Research and Oceanographic Abstracts 15, no. 4
       (August 1968): 497-503. doi:10.1016/0011-7471(68)90057-0.
    '''
    T = (A + B*log(R) + C*(log(R))**3)**-1
    return T


def Steinhart_Hart_resistance(T, A, B, C):
    r'''Calculates the resistance of a Thermistor, using its absolute
    temperature and Steinhard-Hart coefficients according to [1]_. Coefficients
    are calibrated for every single thermistor separately, due to large
    variations in production methods.

    .. math::
        \frac{1}{T} = A + B(\ln R) + C(\ln R)^3

    Parameters
    ----------
    T : float
        Temperature, [K]
    A, B, C : floats
        Coefficients for Steinhart-Hart function, [-]

    Returns
    -------
    R : float
        Resistance of thermistor at T [ohm]

    Notes
    -----
    No warnings produced if outside the suggested range.
    Thermistors are normally valid only for short ranges under 200 K.
    This equation was produced by Sympy, and matches reverse calculations.
    Examples
    --------
    Coefficients for Omega 44004 Thermistor, 2252 ohm at STP:

    >>> Steinhart_Hart_resistance(2252, A=1.468E-3, B=2.383E-4, C=1.007E-7)
    2253.033420231679

    References
    ----------
    .. [1] Steinhart, John S., and Stanley R. Hart. "Calibration Curves for
       Thermistors." Deep Sea Research and Oceanographic Abstracts 15, no. 4
       (August 1968): 497-503. doi:10.1016/0011-7471(68)90057-0.
    '''
    R  = exp(B/(3*C*(sqrt(B**3/(27*C**3) + (A*T - 1)**2/(4*C**2*T**2)) +
    (A*T - 1)/(2*C*T))**(1/3)) - (sqrt(B**3/(27*C**3) +
    (A*T - 1)**2/(4*C**2*T**2)) + (A*T - 1)/(2*C*T))**(1/3))
    return R


def beta_temperature(R, R0, B=3750., T0=298.15):
    r'''Calculates the temperature of a thermistor, using its absolute
    resistance and beta coefficient according to [1]_. Coefficients
    are calibrated for every single thermistor separately, due to large
    variations in production methods. Default value is provided as an example
    only. Also necessary is a reference temperature, and resistivity at the
    reference temperature.

    .. math::
        T = \frac{B}{\log[\frac{R_T}{R_{T0}}\exp(B/T_0)]}

    Parameters
    ----------
    R : float
        Resistance of thermistor at T [ohm]
    R0 : float
        Resistance of thermistor at T0 [ohm]
    B : float, optional
        Beta coefficient, [K]
    T0 : float, optional
        Reference temperature, [K]

    Returns
    -------
    T : float
        Temperature, [K]

    Notes
    -----
    No warnings produced if outside the suggested range.
    Thermistors are normally valid only for short ranges under 200 K.

    Examples
    --------
    >>> beta_temperature(2.5E4, R0=1E4, B=3799.42)
    278.1500055180414

    References
    ----------
    .. [1] Liptak, Bela G. Instrument Engineers' Handbook, 4E, Volume One:
       Process Measurement and Analysis. CRC Press, 2003.
    '''
    T = B*(log(R/R0*exp(B/T0)))**-1
    return T

def beta_resistance(T, R0, B=3750., T0=298.15):
    r'''Calculates the resistance of a thermistor, using its absolute
    temperature and beta coefficient according to [1]_. Coefficients
    are calibrated for every single thermistor separately, due to large
    variations in production methods. Default value is provided as an example
    only. Also necessary is a reference temperature, and resistivity at the
    reference temperature.

    .. math::
        R_T =R_{TO} \exp\left(B\left(\frac{1}{T} - \frac{1}{T_0}\right)\right)

    Parameters
    ----------
    T : float
        Temperature, [K]
    R0 : float
        Resistance of thermistor at T0 [ohm]
    B : float, optional
        Beta coefficient, [K]
    T0 : float, optional
        Reference temperature, [K]

    Returns
    -------
    R : float
        Resistance of thermistor at T [ohm]

    Notes
    -----
    No warnings produced if outside the suggested range.
    Thermistors are normally valid only for short ranges under 200 K.

    Examples
    --------
    >>> beta_resistance(5+273.15, R0=1E4, B=3799.42, T0=298.15)
    25000.00677460831

    References
    ----------
    .. [1] Liptak, Bela G. Instrument Engineers' Handbook, 4E, Volume One:
       Process Measurement and Analysis. CRC Press, 2003.
    '''
    R = R0*exp(B*(1./T - 1./T0))
    return R
