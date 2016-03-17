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



def IAE(ys, ysp):
    r'''Calculates integral absolute error (IAE) for a variable. Also
    called the Sum of Absolute Errors (SAE) .

    .. math::
        IAE = \int_0^\infty |y_{sp}(t) - y_s(t)|dt

    Parameters
    ----------
    ys : array-like
        Varying values of y
    ysp : float
        Desired value for y ("set point")

    Returns
    -------
    IAE : float
        Integral squared error

    Examples
    --------
    >>> IAE([9.9, 9.92, 9.94, 9.96], 10)
    0.27999999999999936
    '''
    iae = sum([abs(yi-ysp) for yi in ys])
    return iae


def ISE(ys, ysp):
    r'''Calculates integral squared error (ISE) for a variable. Also
    called the Sum of Squared Errors (SSE) .

    .. math::
        ISE = \int_0^\infty |y_{sp}(t) - y_s(t)|^2dt

    Parameters
    ----------
    ys : array-like
        Varying values of y
    ysp : float
        Desired value for y ("set point")

    Returns
    -------
    IAE : float
        Integral squared error

    Examples
    --------
    >>> IAE([9.9, 9.92, 9.94, 9.96], 10)
    0.021599999999999935
    '''
    ise = sum([abs(yi-ysp)**2 for yi in ys])
    return ise


def ITAE(ys, ts, ysp):
    r'''Calculates Integral Time Absolute Error (ITAE) for a variable.

    .. math::
        ITAE = \int_0^\infty t|y_{sp}(t) - y_s(t)|dt

    Parameters
    ----------
    ys : array-like
        Varying values of y
    ts : array-like
        Varying values of time, [time]
    ysp : float
        Desired value for y ("set point")

    Returns
    -------
    ITAE : float
        Integral absolute error

    Examples
    --------
    >>> ITAE([9.9, 9.92, 9.94, 9.96], [0, 1, 2, 3], 10)
    0.3199999999999985
    '''
    itae = sum([ts[i]*abs(ys[i]-ysp) for i in range(len(ys))])
    return itae


def ITSE(ys, ts, ysp):
    r'''Calculates Integral Time Squared Error (ITSE) for a variable.

    .. math::
        ITSE = \int_0^\infty t|y_{sp}(t) - y_s(t)|^2 dt

    Parameters
    ----------
    ys : array-like
        Varying values of y
    ts : array-like
        Varying values of time, [time]
    ysp : float
        Desired value for y ("set point")

    Returns
    -------
    ITSE : float
        Integral squared error

    Examples
    --------
    >>> ITAE([9.9, 9.92, 9.94, 9.96], [0, 1, 2, 3], 10)
    0.018399999999999927
    '''
    itse = sum([ts[i]*abs(ys[i]-ysp)**2 for i in range(len(ys))])
    return itse

### Controller form conversions

def interactive_to_noninteractive(Kc, tauI, tauD):
    r'''Converts `Kc`, `tauI`, and `tauD` from the interactive (older)
    controller form to the ideal or standard form (newer). Primes designate
    interactive terms.

    .. math::
        K_c = K_c'\left(\frac{\tau_I' + \tau_D'}{\tau_I'}\right)

        \tau_I = \tau_I' + \tau_D'

        \tau_D = \frac{\tau_I' \tau_D'}{\tau_i' + \tau_D'}

    Parameters
    ----------
    Kc : float
        Controller gain, []
    tauI : float
        Integral time constant, [time]
    tauD : float
        Derivative time constant, [time]

    Returns
    -------
    Kc : float
        Controller gain, []
    tauI : float
        Integral time constant, [time]
    tauD : float
        Derivative time constant, [time]

    Examples
    --------
    >>> interactive_to_noninteractive(0.07236, 7.236, 2.764)
    (0.09999999999999999, 10.0, 2.0000303999999995)
    '''
    Kc_non = Kc*(tauI + tauD)/tauI
    tauI_non = tauI + tauD
    tauD_non = (tauI*tauD)/(tauI + tauD)
    return Kc_non, tauI_non, tauD_non

#print interactive_to_noninteractive(0.07236, 7.236, 2.764)

def noninteractive_to_interactive(Kc, tauI, tauD):
    r'''Converts `Kc`, `tauI`, and `tauD` from the ideal or standard controller
    form (newer) to the interactive (older) form. Primes designate interactive
    terms.

    .. math::
        K_c' = K_c\left(0.5 + \sqrt{0.25 - \tau_D/\tau_I}\right)

        \tau_I' = \tau_I\left(0.5 + \sqrt{0.25 - \tau_D/\tau_I}\right)

        \tau_D' = \frac{\tau_D}{0.5 + \sqrt{0.25 - \tau_D/\tau_I}}

    Parameters
    ----------
    Kc : float
        Controller gain, []
    tauI : float
        Integral time constant, [time]
    tauD : float
        Derivative time constant, [time]

    Returns
    -------
    Kc : float
        Controller gain, []
    tauI : float
        Integral time constant, [time]
    tauD : float
        Derivative time constant, [time]

    Notes
    -----
    If tauD/tauI < 0.25, the conversion is not possible - the square root of a
    negative number is being taken.

    Examples
    --------
    >>> noninteractive_to_interactive(0.1, 10., 2.)
    (0.0723606797749979, 7.23606797749979, 2.7639320225002106)
    '''
    if tauD/tauI > 0.25:
        raise Exception(ValueError)
    Kc_int = Kc*(0.5 + (0.25 - tauD/tauI)**0.5)
    tauI_int = tauI*(0.5 + (0.25 - tauD/tauI)**0.5)
    tauD_int = tauD/(0.5 + (0.25 - tauD/tauI)**0.5)
    return Kc_int, tauI_int, tauD_int


def noninteractive_to_parallel(Kc, tauI, tauD):
    r'''Converts `Kc`, `tauI`, and `tauD` from the ideal or standard
    controller form to `Kp`, `KI`, and `KD` in the parallel form.

    .. math::
        K_p = K_c

        K_I = \frac{K_c}{\tau_I}

        K_D = K_c\tau_D

    Parameters
    ----------
    Kc : float
        Controller gain, []
    tauI : float
        Integral time constant, [time]
    tauD : float
        Derivative time constant, [time]

    Returns
    -------
    Kp : float
        Proportional gain, []
    KI : float
        Integral time constant, [1/time]
    KD : float
        Derivative time constant, [time]

    Examples
    --------
    >>> noninteractive_to_parallel(0.1, 10., 2.)
    (0.1, 0.01, 0.2)
    '''
    Kp = Kc
    KI = Kc/tauI
    KD = Kc*tauD
    return Kp, KI, KD

#print noninteractive_to_parallel(0.1, 10., 2.)

def parallel_to_noninteractive(Kp, KI, KD):
    r'''Converts `Kp`, `KI`, and `KD` in the parallel controller form to
    `Kc`, `tauI`, and `tauD` in the ideal or standard controller form.

    .. math::
        K_c = K_p

        \tau_I = K_c/K_I

        \tau_D = K_D/K_c

    Parameters
    ----------
    Kp : float
        Proportional gain, []
    KI : float
        Integral time constant, [1/time]
    KD : float
        Derivative time constant, [time]

    Returns
    -------
    Kc : float
        Controller gain, []
    tauI : float
        Integral time constant, [time]
    tauD : float
        Derivative time constant, [time]

    Examples
    --------
    >>> parallel_to_noninteractive(0.1, 0.01, 0.2)
    (0.1, 10.0, 2.0)
    '''
    Kc = Kp
    tauI = Kc/KI
    tauD = KD/Kc
    return Kc, tauI, tauD

#print parallel_to_noninteractive(0.1, 0.01, 0.2)

def parallel_to_interactive(Kp, KI, KD):
    r'''Converts Kp`, `KI`, and `KD` in the parallel controller form to `Kc`,
    `tauI`, and `tauD` from the interactive (older) controller form.

    Parameters
    ----------
    Kp : float
        Proportional gain, []
    KI : float
        Integral time constant, [1/time]
    KD : float
        Derivative time constant, [time]

    Returns
    -------
    Kc : float
        Controller gain, []
    tauI : float
        Integral time constant, [time]
    tauD : float
        Derivative time constant, [time]

    Notes
    -----
    Utilizes the functions parallel_to_noninteractive and
    noninteractive_to_interactive to perform the conversions.
    If tauD/tauI < 0.25, in the intermediary step, the conversion is not
    possible - the square root of a negative number is being taken.

    Examples
    --------
    >>> parallel_to_interactive(0.1, 0.01, 0.2)
    (0.0723606797749979, 7.23606797749979, 2.7639320225002106)
    '''
    Kc, tauI, tauD = parallel_to_noninteractive(Kp, KI, KD)
    Kc, tauI, tauD = noninteractive_to_interactive(Kc, tauI, tauD)
    return Kc, tauI, tauD

#print parallel_to_interactive(0.1, 0.01, 0.2)


def interactive_to_parallel(Kc, tauI, tauD):
    r'''Converts `Kc`, `tauI`, and `tauD` from the interactive (older)
    controller form to `Kp`, `KI`, and `KD` in the parallel form.

    Parameters
    ----------
    Kc : float
        Controller gain, []
    tauI : float
        Integral time constant, [time]
    tauD : float
        Derivative time constant, [time]

    Returns
    -------
    Kp : float
        Proportional gain, []
    KI : float
        Integral time constant, [1/time]
    KD : float
        Derivative time constant, [time]

    Notes
    -----
    Utilizes the functions interactive_to_noninteractive and
    noninteractive_to_parallel to perform the conversions.
    If tauD/tauI < 0.25, in the intermediary step, the conversion is not
    possible - the square root of a negative number is being taken.

    Examples
    --------
    >>> interactive_to_parallel(0.07236, 7.236, 2.764)
    (0.09999999999999999, 0.009999999999999998, 0.20000303999999994)
    '''
    Kc, tauI, tauD = interactive_to_noninteractive(Kc, tauI, tauD)
    Kp, KI, KD = noninteractive_to_parallel(Kc, tauI, tauD)
    return Kp, KI, KD

#print interactive_to_parallel(0.07236, 7.236, 2.764)


### Controller form time unit conversions