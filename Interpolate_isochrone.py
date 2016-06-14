__author__ = 'Felipe'

from Plot_isochrones_v2 import *
import numpy as np
import scipy.interpolate as intpl
from matplotlib import pyplot as plt

def interp_isoc(Te, Mv, m, Te_lim = 0, Mv_lim = 0, nintpl = 1000):
    if Te_lim and Mv_lim:
        filter = (Te >= Te_lim[0]) & (Te < Te_lim[1]) & (Mv >= Mv_lim[0]) & (Mv < Mv_lim[1])
        Te = Te[filter]
        Mv = Mv[filter]
        m = m[filter]

    Te_fun = intpl.interp1d(m, Te)
    Mv_fun = intpl.interp1d(m, Mv)

    m_intpl = np.linspace(m[0], m[-1], nintpl)
    Te_intpl = Te_fun(m_intpl)
    Mv_intpl = Mv_fun(m_intpl)

    return(Te_intpl, Mv_intpl, m_intpl)

def interp_isoc_wstage(Te, Mv, m, stage, nintpl = 1000):
    Te_intpl = np.array(())
    Mv_intpl = np.array(())
    m_intpl = np.array(())

    for s in range(9):
        if s in stage:
            filter = stage == s
            if filter.sum() > 2:
                Te_stage = Te[filter]
                Mv_stage = Mv[filter]
                m_stage = m[filter]

                Te_temp, Mv_temp, m_temp = interp_isoc(Te_stage, Mv_stage, m_stage, nintpl = nintpl)

            else:
                Te_temp = Te[filter]
                Mv_temp = Mv[filter]
                m_temp = m[filter]

            Te_intpl = np.concatenate((Te_intpl, Te_temp))
            Mv_intpl = np.concatenate((Mv_intpl, Mv_temp))

    return Te_intpl, Mv_intpl, m_intpl

if 0:
    data = np.loadtxt('C:\Users\Felipe\Documents\Doutorado\Projeto Idade Isocronal\Isochrone\Isochrones_set_parsec1.0_res0.5Ga\parsec1.0_res0.5_3.00e9.dat', usecols=(0,2,5,7,17))
    Z = data[:,0]
    m = data[:,1]
    Te = 10**data[:,2]
    Mv = data[:,3]
    s = data[:,4]

    Z0 = 0.0025
    filter = Z == Z0

    m = m[filter]
    Te = Te[filter]
    Mv = Mv[filter]
    s = s[filter]

    plt.plot(Te, Mv, c = 'b')
#    plt.scatter(Te, Mv, c = 'b')

    Te_intpl, Mv_intpl, m_intpl = interp_isoc(Te, Mv, m, nintpl = 601)
#    plt.scatter(Te_intpl, Mv_intpl, c = 'r')
    plt.plot(Te_intpl, Mv_intpl, '--r')

    Te_intpl, Mv_intpl, m_intpl = interp_isoc_wstage(Te, Mv, m, stage = s, nintpl = 500)
#    plt.scatter(Te_intpl, Mv_intpl, c = 'g')
    plt.plot(Te_intpl, Mv_intpl, '--g')

    ax = plt.subplot(111)
    ax.invert_xaxis()
    ax.invert_yaxis()

    plt.show()

