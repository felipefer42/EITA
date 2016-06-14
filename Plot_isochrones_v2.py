__author__ = 'Felipe'

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib import gridspec

def get_filename(t, name_base = 'set1_res0.01', name_type=1):
    if name_type == 1:
        t_str = "_%0.2f"%(t/1e9)+'e9'
    elif name_type == 2:
        t_str = "_t%0.1f"%(t/1e9)
    return(name_base+t_str+'.dat')

def load_isochrone_ZmTeMv(t, path, name_base = 'set1_res0.01'):
    filename = get_filename(t, name_base)
    data = np.loadtxt(path+'\\'+filename, usecols=(0,2,5,7))
    Z = data[:,0]
    m = data[:,1]
    Te = 10**data[:,2]
    Mv = data[:,3]

    return(Z, m, Te, Mv)

def config_cbar(t, cmap = 'cool'):
    cmx = plt.get_cmap(cmap)
    X = [[0,0], [0,0]]
    levels = np.arange(t.min(), t.max(), t[1]-t[0])
    CS3 = plt.contourf(X, levels, cmap = cmx)

    plt.clf()
    return CS3

def plot_isochrone(t, path, name_base, Te_lim = (4000., 8000.), Mv_lim = (2., 5.5), Z_isoc = 0.002, cmap = 'cool', alpha = 1, invert_axis = True, rep_mass = True, ax = plt):

    if ax == plt:
        ax = plt.subplot(111)

    # Configuring colormap
    cmx = plt.get_cmap(cmap)
    cNorm = colors.Normalize(vmin = t.min(), vmax = t.max())
    scalarMap = cm.ScalarMappable(norm = cNorm, cmap = cmx)

    # Ploting isochrone relative to each age in t
    for i in range(len(t)):
        Z, m, Te, Mv = load_isochrone_ZmTeMv(t[i], path, name_base)
        filter = Z == Z_isoc
        Te = Te[filter]
        Mv = Mv[filter]

        color = scalarMap.to_rgba(t[i])
        ax.plot(Te, Mv, color = color, alpha = alpha)

        if rep_mass == 'no_dots':
            pass
        else:
            if rep_mass:
                m = m[filter]
                ax.scatter(Te, Mv, s = 0.2*np.exp(5*m), c = color, alpha = alpha)
            else:
                ax.scatter(Te, Mv, s = 5, c = color, edgecolor = color, alpha = alpha)

    # Setting limits
    ax.set_xlim(Te_lim)
    ax.set_ylim(Mv_lim)

    # Inverting axis
    if invert_axis:
        ax.invert_xaxis()
        ax.invert_yaxis()

Teste_plot_isochrone = False
if Teste_plot_isochrone:
    t = np.arange(1e9, 13e9, 1e9)
    path = r'C:\Users\Felipe\Documents\Doutorado\Projeto Idade Isocronal\Isochrone\Isochrones_setMarigo08_res0.5Ga'
    name_base = r'Marigo08_res0.5'

    config_cbar(t)
    plot_isochrone(t, path, name_base, Te_lim = (5000,9000), Mv_lim = (1, 5))

    plt.show()

def plot_n_ischrones(t, path, name_base, Te_lim = (4000., 8000.), Mv_lim = (2., 5.5), Z_isocs = (0.0005, 0.001, 0.002, 0.003, 0.0035), cmap = 'cool', alpha = (1, 1, 1, 1, 1),invert_axis = True, rep_mass = True, figsize = (16,4), cbar_ratio = 15):

    n = len(Z_isocs)

    CS3 = config_cbar(t)

#    f, ax = plt.subplots(1, (n+1), sharey=False, figsize = (16,4))
    f = plt.figure(figsize = figsize)
    gs = gridspec.GridSpec(1, (n+1), width_ratios=(n*[cbar_ratio]+[1]))
    ax = (n+1)*[0]
    for i in range(n+1):
        ax[i] = plt.subplot(gs[i])

    for i in range(n):
        plot_isochrone(t = t, path = path, name_base = name_base, Te_lim = Te_lim, Mv_lim = Mv_lim, Z_isoc = Z_isocs[i], cmap = cmap, alpha = alpha[i], invert_axis = invert_axis, rep_mass = rep_mass, ax = ax[i])

    cbar = plt.colorbar(CS3, cax = ax[-1])

    for i in range(n):
        ax[i].set_xlabel('Te (K)')
    ax[0].set_ylabel('Mv')

    for i in range(1,n):
        ax[i].set_yticklabels([])

    return(f, ax, cbar)

##########################################################
# Teste

if 0:
    t = np.arange(0.5e9, 13.51e9, 0.5e9)
    path = r'C:\Users\Felipe\Documents\Doutorado\Projeto Idade Isocronal\Isochrone\Isochrones_set_parsec1.0_res0.5Ga'
    name_base = r'parsec1.0_res0.5'

    Te_obs = 5777
    Te_err = 100
    Mv_obs = 4.75
    Mv_err = 0.1
    Z_obs = 0.02
    Z_err = 0.005

    Z_isocs = [Z_obs-2*Z_err, Z_obs-Z_err, Z_obs, Z_obs + Z_err, Z_obs+2*Z_err]
    Z_text = [r'$Z = {0} - 2\,\sigma_Z$'.format(Z_isocs[2]),
              r'$Z = {0} - \sigma_Z$'.format(Z_isocs[2]),
              r'$Z = {0}$'.format(Z_isocs[2]),
              r'$Z = {0} + \sigma_Z$'.format(Z_isocs[2]),
              r'$Z = {0} + 2\,\sigma_Z$'.format(Z_isocs[2])]
    Alpha = [1, 1, 1, 1, 1]

    f, ax, cbar = plot_n_ischrones(t, path, name_base, Te_lim = (5000, 6500), Mv_lim = (3.5, 5.5), Z_isocs = Z_isocs, rep_mass= False)
    for i in range(5):
        ax[i].xaxis.set_ticks(np.arange(5000, 6500, 250))

    for i in range(5):
        ax[i].errorbar(Te_obs, Mv_obs, xerr = Te_err, yerr = Mv_err, fmt = 'o', ecolor = 'k')
        ax[i].text(6250, 5.25, Z_text[i])

    cbar.set_label('age')

    # Ploting points
    plt.tight_layout()

    plt.subplots_adjust(wspace = 0, hspace = 0)
    plt.show()

# Single plot
if 0:
    t = np.arange(0.5e9, 13.51e9, 0.5e9)
    path = r'C:\Users\Felipe\Documents\Doutorado\Projeto Idade Isocronal\Isochrone\Isochrones_set_parsec1.0_res0.1Ga'
    name_base = r'parsec1.0_res0.5'

    Z_obs = 0.0075

    Z_isocs = [Z_obs]
    Z_text = [r'$Z = {0}$'.format(Z_obs)]
    Alpha = [1]
    figsize = (9,7.5)

    f, ax, cbar = plot_n_ischrones(t, path, name_base, Te_lim = (4800, 6300), Mv_lim = (3.5, 6.5), Z_isocs = Z_isocs, rep_mass= False, figsize = figsize, cbar_ratio = 20)
    for i in range(1):
        ax[i].xaxis.set_ticks(np.arange(5000, 6300, 500))

#    for i in range(5):
#        ax[i].errorbar(Te_obs, Mv_obs, xerr = Te_err, yerr = Mv_err, fmt = 'o', ecolor = 'k')
#        ax[i].text(6250, 5.25, Z_text[i])

    cbar.set_label('age')

    # Ploting points
    plt.tight_layout()

    plt.subplots_adjust(wspace = 0, hspace = 0)

    # Takeda+2007
    obj_id = ['HD101614', 'HD193307', 'HD120237', 'HD10700', 'HD144628']
    Te = np.array([5687.07, 6033.01, 6126.96, 5282.91, 5021.01])
    Mv = np.array([4.31, 3.73, 4.35, 5.58, 6.1])

    #Z_obs = np.array([0.027, 0.009, 0.014, 0.006, 0.007])

    ax[0].text(6000, 6, "Z_isoc = {0}".format(Z_obs))

    Mv_err = [0.1, 0.2, 0.3]
    Te_err = [50, 100, 200]

    for i in range(len(obj_id)):
        ax[0].text(Te[i]+15, Mv[i]-0.05, obj_id[i])

    ax[0].scatter(Te, Mv, c = ['b', 'b', 'b', 'b', 'r'])

    #ax[0].errorbar(Te, Mv, xerr = 0, yerr = 0, fmt = 'o', ecolor = 'k')
    #ax[0].errorbar([7400, 7050, 6500], [5,5,5], xerr = Te_err, yerr = Mv_err, fmt = '.', ecolor = 'k')

    #ax[0].text(7400-20, 5-0.1, r'$1\sigma$')
    #ax[0].text(7050-20, 5-0.1, r'$2\sigma$')
    #ax[0].text(6500-20, 5-0.1, r'$4\sigma$')

    plt.show()

if 0:
    t = np.arange(0.5e9, 13.51e9, 0.5e9)
    path = r'C:\Users\Felipe\Documents\Doutorado\Projeto Idade Isocronal\Isochrone\Isochrones_set_parsec1.0_res0.5Ga'
    name_base = r'parsec1.0_res0.5'

    Z_obs = 0.02

    Z_isocs = [Z_obs]
    Z_text = [r'$Z = {0}$'.format(Z_obs)]
    Alpha = [1]
    figsize = (7.612,6.792)

    f, ax, cbar = plot_n_ischrones(t, path, name_base, Te_lim = (4400, 7600), Mv_lim = (0.9, 6.1), Z_isocs = Z_isocs, rep_mass= 'no_dots', figsize = figsize, cbar_ratio = 20)
    for i in range(1):
        ax[i].xaxis.set_ticks(np.arange(4500, 8000, 500))

#    for i in range(5):
#        ax[i].errorbar(Te_obs, Mv_obs, xerr = Te_err, yerr = Mv_err, fmt = 'o', ecolor = 'k')
#        ax[i].text(6250, 5.25, Z_text[i])

    cbar.set_label('')

    # Ploting points
    plt.tight_layout()

    plt.subplots_adjust(wspace = 0, hspace = 0)
    # Test stars:
    Mv = [6, 5.5, 5, 5, 5, 4.5, 4.5, 4.5, 4, 4, 4, 4, 3.5, 3.5, 3.5, 3.5, 2.5, 2.5, 2.5, 2.5, 2.5, 1.5, 1.5, 1.5, 1.5, 1.5]
    Te = [5000,5300, 5450, 5550, 5650, 5600, 5750, 5900, 5000, 5420, 5840, 6250, 4850, 5430, 6010, 6600, 4800, 5435, 6070, 6705, 7350, 4700, 5375, 6050, 6725, 7400]
    obj_id = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    obj_id.reverse()

    Mv_err = [0.1, 0.2, 0.3]
    Te_err = [50, 100, 200]

    #for i in range(len(obj_id)):
    #    ax[0].text(Te[i]+15, Mv[i]-0.05, obj_id[i])

    #ax[0].errorbar(Te, Mv, xerr = 0, yerr = 0, fmt = 'o', ecolor = 'k')
    #ax[0].errorbar([4750], [5], xerr = Te_err[0], yerr = Mv_err[0], fmt = '.', ecolor = 'k')

    #ax[0].text(4750-20, 5-0.1, r'$1\sigma$')

    plt.savefig('C:\Users\Felipe\Documents\Doutorado\Apresentacoes\Isocronas_exemplo2.png', dpi = 300)
    plt.show()
