# -*- coding: utf-8 -*-

import clamped_spline as csp
import math

def get_curve(plazos, tasas):
    """
    Dados numpy arrays plazos (en años) y tasas se retorna
    un clamped cubic spline que permite interpolar en la
    curva a cualquier plazo y también obtener
    las derivadas primera y segunda en ese plazo.
    """
    return csp.get_clamped_spline_interpol(plazos, tasas)


def fwd(zrate, t):
    """
    Tasa forward instantánea en t
    t: tiempo en años
    """
    return zrate(t)[0] + t * zrate(t)[1]


def dfwd(zrate, t):
    """
    Derivada de la tasa forward instantánea
    t: tiempo en años
    """
    return 2 * zrate(t)[1] + t * zrate(t)[2]


def theta(zrate, t, sigma, gamma):
    aux = (sigma ** 2) / (2.0 * gamma) * (1 - math.exp(-2.0 * gamma * t))
    return dfwd(zrate, t) + gamma * fwd(zrate, t) + aux


def b_hw(gamma, t, T):
    """
    Coeficiente B del modelo de HW
    gamma : parámetro de reversión 
    t : plazo de....
    T: plazo ....
    """
    aux = 1 - math.exp(- gamma * (T - t))
    return aux / gamma


def a_hw(zrate, gamma, sigma, t, T, verbose = False):
    """
    Coeficiente A del modelo de HW.
    verbose: cuando es True imprime los valores de c1, c2 y c3.
    """
    b = b_hw(gamma, t, T)
    dfT = math.exp(-zrate(T)[0] * T)
    dft = math.exp(-zrate(t)[0] * t)
    c1 = math.log(dfT / dft)
    c2 = b * fwd(zrate, t)
    c3 = (sigma**2) / (4 * gamma) * (b**2) * (1 - math.exp(-2 * gamma * t))
    if verbose:
        print("c1: " + str(c1))
        print("c2: " + str(c2))
        print("c3: " + str(c3))
    return c1 + c2 - c3
