"""
Copyright (C) 2019-2020 Emanuele Paci, Simon P. Skinner, Michele Stofella

This program is free software: you can redistribute it and/or modify
it under the terms of version 2 of the GNU General Public License as published
by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from math import exp, log10

# Parameters for a protonated protein in D2O
# Here are the parameters measured in Bai et al. (1993)
# The parameters for D and E are based on the work of Mori and al. (1997)
# and measured by and calculated in the functions acid and base
para = {
    "A": [0.00,   0.00,  0.00,  0.00],
    "C": [-0.54, -0.46,  0.62,  0.55],
    "F": [-0.52, -0.43, -0.24,  0.06],
    "G": [-0.22,  0.22, -0.03,  0.17],
    "I": [-0.91, -0.59, -0.73, -0.23],
    "K": [-0.56, -0.29, -0.04,  0.12],
    "L": [-0.57, -0.13, -0.58, -0.21],
    "M": [-0.64, -0.28, -0.01,  0.11],
    "N": [-0.58, -0.13,  0.49,  0.32],
    "P": [99999, -0.19, 99999, -0.24], # trans
    "B": [99999, -0.85, 99999,  0.60], # proline cis
    "Q": [-0.47, -0.27,  0.06,  0.20],
    "R": [-0.59, -0.32,  0.08,  0.22],
    "S": [-0.44, -0.39,  0.37,  0.30],
    "T": [-0.79, -0.47, -0.07,  0.20],
    "V": [-0.74, -0.30, -0.70, -0.14],
    "W": [-0.40, -0.44, -0.41, -0.11],
    "Y": [-0.41, -0.37, -0.27,  0.05],
}

rho_Nterm_acid = -1.32
rho_Nterm_base = 1.62

lamb_Cterm_acid = 0.96
lamb_Cterm_base = -1.80

pKD = 15.05
R = 1.987

ka_pdla = 10**(1.62) / 60
kb_pdla = 10**(10.18) / 60
kw_pdla = 10**(-1.5) / 60

ka = 10**(2.04) / 60
kb = 10**(10.36) / 60
kw = 10**(-1.5) / 60

Ea = 14000
Eb = 17000
Ew = 19000


def get_D(pH):
    return 10**(-pH)


def get_OD(pH):
    return 10**(pH - pKD)


def get_temperature_normalization(temperature):
    return (1 / temperature - 1 / 293.15) / R


def get_pK_his(temperature):
    Ea_his = 7500
    return -log10(10**(-7.42) * exp(-Ea_his * (1 / temperature - 1 / 278.15) / R))


def get_pK_asp(temperature):
    Ea_asp = 1000
    return -log10(10**(-4.48) * exp(-Ea_asp * (1 / temperature - 1 / 278.15) / R))


def get_pK_glu(temperature):
    Ea_glu = 1083
    return -log10(10**(-4.93) * exp(-Ea_glu * (1 / temperature - 1 / 278.15) / R))


def get_Fta(temperature):
    return exp(-Ea * get_temperature_normalization(temperature))


def get_Ftb(temperature):
    return exp(-Eb * get_temperature_normalization(temperature))


def get_Ftw(temperature):
    return exp(-Ew * get_temperature_normalization(temperature))
