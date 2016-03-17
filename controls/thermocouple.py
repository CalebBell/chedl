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
from math import exp
from controls.sensor import Sensor
from numpy import roots, real, imag
from scipy.constants import *


def calc_EMF(T, coefs_splits, coefs_piecewise, coefs_exponential):
    t = T - 273.15
    for i in range(len(coefs_splits)-1):
        if T < coefs_splits[i+1]:
            break

    E = sum([c*t**j for j, c in enumerate(coefs_piecewise[i])])
    if coefs_exponential and T > 273.15:
        # Defined for K thermocouple only, otherwise a generic definition
        # would be needed
        E += coefs_exponential[0]*exp(coefs_exponential[1]*((t-coefs_exponential[2])**2))
    return E

def calc_T(E, Tref, coefs_splits, coefs_piecewise, coefs_exponential):
    Eref = calc_EMF(Tref, coefs_splits, coefs_piecewise, coefs_exponential)
    # Brute force the polynomials to determine the solution
    for i in range(len(coefs_splits)-1):
        coefs = coefs_piecewise[i]
        Tmin_coefs, Tmax_coefs = coefs_splits[i], coefs_splits[i+1]

        coefs.reverse()
        coefs[-1] -= E + Eref
        ans = roots(coefs)
        for j in ans:
            if not imag(j) and real(j) + 273.15 > Tmin_coefs and real(j) + 273.15 < Tmax_coefs:
                return real(j) + 273.15
    # If the above didn't work, T is outside of safe limits.
    # TODO: Try the lower limit, then the higher limit, with fsolve.
    # No guarantees though.
    return None

#def sum_emfs(t, coefs):
#    E = sum([c*t**i for i, c in enumerate(coefs)])
#    return E

#def thermocouple(V=None, T=None, Tref=298.15, kind='K'):
##    T = T - 273.15
#    if kind == 'K':
#        ais = [0.118597600000E+00, -0.118343200000E-03, 0.126968600000E+03]
#        if T < 0.000:
#            if T < -270:
#            coefs = [0.000000000000E+00, 0.394501280250E-01,
#                     0.236223735980E-04, -0.328589067840E-06,
#                     -0.499048287770E-08, -0.675090591730E-10,
#                     -0.574103274280E-12, -0.310888728940E-14,
#                     -0.104516093650E-16, -0.198892668780E-19,
#                     -0.163226974860E-22]
#            E = sum_emfs(T, Tref, coefs)
#        else:
#            coefs = [-0.176004136860E-01, 0.389212049750E-01,
#                     0.185587700320E-04, -0.994575928740E-07,
#                     0.318409457190E-09, -0.560728448890E-12,
#                     0.560750590590E-15, -0.320207200030E-18,
#                     0.971511471520E-22, -0.121047212750E-25]
#            E = sum_emfs(T, Tref, coefs) + ais[0]*exp(ais[1]*((T-ais[2])**2)) - ais[0]*exp(ais[1]*((Tref-ais[2])**2))
#    return E

#print thermocouple(T=70.6, Tref=21.3)


Thermocouple_list = ['Au-Pt',  'E',  'G',  'IrRh 0.4 vs Ir',  'J',  'K',  'K vs AuFe 0.07',  'M',  'N',  'P',  'PtMo 0.05 vs PtMo 0.001',  'Pt-Pd',  'PtRh 0.4 vs PtRh 0.2',  'R',  'S',  'T',  'B']

class Thermocouple(Sensor):
    '''RTD Class'''
    def __init__(self,
                 variable_units='K', signal_units='mV',
                 deadtime=0, offset=0, output_rounding=None,
                 method='K', Tref=298.15):
        self.method = method
        self.Tref = Tref
        self.signal_to_variable_function = calc_T
        self.variable_to_signal_function = calc_EMF
        self.coefs_exponential = None

        if self.method == 'B':
            self.coefs_piecewise = [[0.000000000000E+00, -0.246508183460E-03, 0.590404211710E-05, -0.132579316360E-08, 0.156682919010E-11, -0.169445292400E-14, 0.629903470940E-18],
                                    [-0.389381686210E+01, 0.285717474700E-01, -0.848851047850E-04, 0.157852801640E-06, -0.168353448640E-09, 0.111097940130E-12, -0.445154310330E-16, 0.989756408210E-20, -0.937913302890E-24]]
            self.coefs_splits = [273.15, 630.615 + 273.15, 1820.000 + 273.15]
            self.description = 'Pt-30% Rh versus Pt-6% Rh'
        elif self.method == 'E':
            self.coefs_piecewise = [[0.000000000000E+00, 0.586655087080E-01, 0.454109771240E-04, -0.779980486860E-06, -0.258001608430E-07, -0.594525830570E-09, -0.932140586670E-11, -0.102876055340E-12, -0.803701236210E-15, -0.439794973910E-17, -0.164147763550E-19, -0.396736195160E-22, -0.558273287210E-25, -0.346578420130E-28],
                                    [0.000000000000E+00, 0.586655087100E-01, 0.450322755820E-04, 0.289084072120E-07, -0.330568966520E-09, 0.650244032700E-12, -0.191974955040E-15, -0.125366004970E-17, 0.214892175690E-20, -0.143880417820E-23, 0.359608994810E-27]]
            self.coefs_splits = [3.15, 273.15, 1000 + 273.15]
            self.description = 'Ni-Cr alloy versus a Cu-Ni alloy'
        elif self.method == 'J':
            self.coefs_piecewise = [[0.000000000000E+00, 0.503811878150E-01, 0.304758369300E-04, -0.856810657200E-07, 0.132281952950E-09, -0.170529583370E-12, 0.209480906970E-15, -0.125383953360E-18, 0.156317256970E-22],
                                    [0.296456256810E+03, -0.149761277860E+01, 0.317871039240E-02, -0.318476867010E-05, 0.157208190040E-08, -0.306913690560E-12]]
            self.coefs_splits = [63.15, 760 + 273.15, 1200.000 + 273.15]
            self.description = 'Fe versus a Cu-Ni alloy'
        elif self.method == 'K':
            self.coefs_piecewise = [[0.000000000000E+00, 0.394501280250E-01, 0.236223735980E-04, -0.328589067840E-06, -0.499048287770E-08, -0.675090591730E-10, -0.574103274280E-12, -0.310888728940E-14, -0.104516093650E-16, -0.198892668780E-19, -0.163226974860E-22],
                                    [-0.176004136860E-01, 0.389212049750E-01, 0.185587700320E-04, -0.994575928740E-07, 0.318409457190E-09, -0.560728448890E-12, 0.560750590590E-15, -0.320207200030E-18, 0.971511471520E-22, -0.121047212750E-25]]
            self.coefs_exponential = [0.118597600000E+00, -0.118343200000E-03, 0.126968600000E+03]
            self.coefs_splits = [3.15, 273.15, 1372.000 + 273.15]
            self.description = 'Ni-Cr alloy versus Ni-Al alloy'
        elif self.method == 'N':
            self.coefs_piecewise = [[0.000000000000E+00, 0.261591059620E-01, 0.109574842280E-04, -0.938411115540E-07, -0.464120397590E-10, -0.263033577160E-11, -0.226534380030E-13, -0.760893007910E-16, -0.934196678350E-19],
                                    [0.000000000000E+00, 0.259293946010E-01, 0.157101418800E-04, 0.438256272370E-07, -0.252611697940E-09, 0.643118193390E-12, -0.100634715190E-14, 0.997453389920E-18, -0.608632456070E-21, 0.208492293390E-24, -0.306821961510E-28]]
            self.coefs_splits = [3.15, 273.15, 1300.000 + 273.15]
            self.description = 'Ni-Cr-Si alloy versus Ni-Si-Mg alloy'
        elif self.method == 'R':
            self.coefs_piecewise = [[0.000000000000E+00, 0.528961729765E-02, 0.139166589782E-04, -0.238855693017E-07, 0.356916001063E-10, -0.462347666298E-13, 0.500777441034E-16, -0.373105886191E-19, 0.157716482367E-22, -0.281038625251E-26],
                                    [0.295157925316E+01, -0.252061251332E-02, 0.159564501865E-04, -0.764085947576E-08, 0.205305291024E-11, -0.293359668173E-15],
                                    [0.152232118209E+03, -0.268819888545E+00, 0.171280280471E-03, -0.345895706453E-07, -0.934633971046E-14]]
            self.coefs_splits = [223.15, 1064.180 + 273.15, 1664.500 + 273.15, 1768.1 + 273.15]
            self.description = 'Pt-13% Rh versus Pt'
        elif self.method == 'S':
            self.coefs_piecewise = [[0.000000000000E+00, 0.540313308631E-02, 0.125934289740E-04, -0.232477968689E-07, 0.322028823036E-10, -0.331465196389E-13, 0.255744251786E-16, -0.125068871393E-19, 0.271443176145E-23],
                                    [0.132900444085E+01, 0.334509311344E-02, 0.654805192818E-05, -0.164856259209E-08, 0.129989605174E-13],
                                    [0.146628232636E+03, -0.258430516752E+00, 0.163693574641E-03, -0.330439046987E-07, -0.943223690612E-14]]
            self.coefs_splits = [223.15, 1064.180+273.15, 1664.5 + 273.15, 1768.1 + 273.15]
            self.description = 'Pt-10% Rh versus Pt'
        elif self.method == 'T':
            self.coefs_piecewise = [[0.000000000000E+00, 0.387481063640E-01, 0.441944343470E-04, 0.118443231050E-06, 0.200329735540E-07, 0.901380195590E-09, 0.226511565930E-10, 0.360711542050E-12, 0.384939398830E-14, 0.282135219250E-16, 0.142515947790E-18, 0.487686622860E-21, 0.107955392700E-23, 0.139450270620E-26, 0.797951539270E-30],
                                    [0.000000000000E+00, 0.387481063640E-01, 0.332922278800E-04, 0.206182434040E-06, -0.218822568460E-08, 0.109968809280E-10, -0.308157587720E-13, 0.454791352900E-16, -0.275129016730E-19]]
            self.coefs_splits = [3.15, 273.15, 400.000 + 273.15]
            self.description = 'Cu versus a Cu-Ni alloy'

        # ASTM E1751 supplementary polynomials
        elif self.method == 'G':
            self.coefs_piecewise = [[0.0, 0.0012792201, 2.1634754e-05, -1.1393234e-08, 4.3850022e-12, -1.7089202e-15],
                                    [-1.1064412, 0.0094962455, -3.6467516e-06, 3.114133e-08, -3.8615222e-11, 2.4455012e-14, -8.9888053e-18, 1.8120237e-21, -1.5534591e-25]]
            self.coefs_splits = [273.15, 630.615 + 273.15, 2315 + 273.15]
            self.description = 'Tungsten versus Tungsten - 26% Rhenium'
        elif self.method == 'P':
            self.coefs_piecewise = [[0.0, 0.029819716, 3.5175152e-05, -3.4878428e-08, 1.4851327e-11, -3.6375467e-15],
                                    [-8.9621838, 0.0853772, -0.00010570233, 1.5424937e-07, -1.2855115e-10, 5.443876e-14, -9.3211269e-18]]
            self.coefs_splits = [273.15, 746.6 + 273.15, 273.15 + 1395]
            self.description = 'Platinel II'
        elif self.method == 'K vs AuFe 0.07':
            self.coefs_piecewise = [[0.0, 0.022272367466, 3.6406179664e-06, -1.5967928202e-07, -4.5260169888e-09, 4.0432555769e-11, 4.9063035769e-12, 1.2272348484e-13, 1.6829773697e-15, 1.4636450149e-17, 8.4287909747e-20, 3.2146639387e-22, 7.8225430483e-25, 1.1010930596e-27, 6.826366158e-31]]
            self.coefs_splits = [0.15, 273.15+7]
            self.description = 'K vs Gold - 0.07% Iron'
        elif self.method == 'PtMo 0.05 vs PtMo 0.001':
            self.coefs_piecewise = [[0.0, 0.010501456, 2.8410937e-05, -4.3368594e-08, 1.058577e-10, -2.384895e-13, 3.3574252e-16, -2.0186476e-19],
                                    [6.8354086, -0.048776479, 0.00024913353, -4.9920472e-07, 6.4615219e-10, -5.3071212e-13, 2.6865173e-16, -7.6717268e-20, 9.4670862e-24]]
            self.coefs_splits = [273.15, 273.15 + 491., 273.15 + 1600]
            self.description = 'Platinum - 5% Molybdenum versus Platinum - 0.1% Molybdenum'
        elif self.method == 'PtRh 0.4 vs PtRh 0.2':
            self.coefs_piecewise = [[0.0, 0.00036246289, 3.936032e-07, 4.2594137e-10, 1.0382985e-12, -1.5406939e-15, 1.0033974e-18, -2.849716e-22],
                                    [-0.91201877, 0.0035246931, -3.9077442e-06, 3.6728697e-09, -1.082471e-12, 1.151628e-16, -1.261964e-20]]
            self.coefs_splits = [273.15, 951.7 + 273.15, 1888 + 273.15]
            self.description = 'Platinum 40% Rhodium versus Platinum - 20% Rhodium'

        elif self.method == 'M':
            self.coefs_piecewise = [[0.0, 0.03690092195, 4.408522682e-05, -3.142898226e-08, -1.02521613e-10, 1.846977453e-13, -9.738054601e-17, -3.3943879e-19],
                                    [-11.45582129, 0.2059913943, -0.0008846963426, 2.650568429e-06, -4.958763813e-09, 6.145877457e-12, -5.041679909e-15, 2.627522669e-18, -7.864442961e-22, 1.027600874e-25]]
            self.coefs_splits = [-50 + 273.15, 370.8 + 273.15, 1410 + 273.15]
            self.description = 'Nickel-18% Molybdenum vs Nickel-0.008% Cobalt'
        elif self.method == 'IrRh 0.4 vs Ir':
            self.coefs_piecewise = [[0.0, 0.0030870016, 6.9649773e-06, -7.8890504e-09, 2.7700591e-12, 2.6762413e-14, -1.041804e-16, 1.5270867e-19, -7.9634082e-23],
                                    [-0.096839082, 0.0036588615, 5.7455189e-06, -6.0547943e-09, 2.7235393e-12, -5.1797037e-16, 3.0821886e-20]]
            self.coefs_splits = [273.15, 630.615 + 273.15, 2110 + 273.15]
            self.description = 'Iridium - 40% Rhodium vs Iridium'
        elif self.method == 'Au-Pt':
            self.coefs_piecewise = [[0.0, 0.00603619861, 1.93672974e-05, -2.22998614e-08, 3.28711859e-11, -4.24206193e-14, 4.56927038e-17, -3.39430259e-20, 1.4298159e-23, -2.51672787e-27]]
            self.coefs_splits = [273.15, 1273.15]
            self.description = 'Gold vs. Platinum'
        elif self.method == 'Pt-Pd':
            # Coefficients originally in microvolts, converted to mV
            self.coefs_piecewise = [[0.0, 0.005296958, 4.610494e-06, -9.602271e-09, 2.992243e-11, -2.012523e-14, -1.268514e-17, 2.257823e-20, -8.510068e-24],
                                    [-0.4977137, 0.010182545, -1.5793515e-05, 3.63617e-08, -2.6901509e-11, 9.5627366e-15, -1.3570737e-18]]
            self.coefs_splits = [273.15, 660.323 + 273.15, 1500 + 273.15]
            self.description = 'Platinum vs. Palladium'
        else:
            raise Exception(ValueError)

        variable_min = self.coefs_splits[0]
        variable_max = self.coefs_splits[-1]

        signal_min = self.variable_to_signal(variable_min, Tref=self.Tref)
        signal_max = self.variable_to_signal(variable_max, Tref=self.Tref)

        Sensor.__init__(self, variable_min=variable_min, variable_max=variable_max,
                 variable_units=variable_units, signal_min=signal_min,
                 signal_max=signal_max, signal_units=signal_units,
                 deadtime=deadtime, offset=offset)


    def signal_to_variable(self, signal, Tref=None):
        if Tref:
            self.Tref = Tref
        T =  self.signal_to_variable_function(signal, self.Tref, self.coefs_splits, self.coefs_piecewise, self.coefs_exponential)
        return T
    def variable_to_signal(self, variable, Tref=None):
        if Tref:
            self.Tref = Tref
#        print variable, self.Tref
#        print self.coefs_splits, self.coefs_piecewise, self.coefs_exponential
        E = (self.variable_to_signal_function(variable, self.coefs_splits, self.coefs_piecewise, self.coefs_exponential)
            - self.variable_to_signal_function(self.Tref, self.coefs_splits, self.coefs_piecewise, self.coefs_exponential))
        return E

    def variable_function(self, variable):
        E = calc_EMF(variable, self.coefs_splits, self.coefs_piecewise, self.coefs_exponential)
        return E




def plot_thermocouples():
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    from matplotlib.backends.backend_pdf import PdfPages
    fontP = FontProperties()
    fontP.set_size('small')
    fig = plt.figure(figsize=(11.69,8.27))
    ax = plt.subplot(111)
    for ttype in Thermocouple_list:
        t = Thermocouple(method=ttype)
        Tmin, Tmax = t.variable_min, t.variable_max
        Es = []
        Ts = np.linspace(Tmin, Tmax, 1000)
        for T in Ts:
            Es.append(t.variable_function(T))
        ax.plot(Ts, Es, label=t.method)
        ax.text(Ts[-1], Es[-1], t.method)
    box = ax.get_position()
#    ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
#    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop=fontP)
    plt.title('EMF vs temperature for most Thermocouple types')
    plt.xlabel('Temperature, K')
    plt.ylabel('EMF in mV')
    plt.grid(True)
    ax.xaxis.grid
    plt.show()
    pp = PdfPages('EMF vs T.pdf')
    plt.savefig(pp, format='pdf')
    pp.close()

