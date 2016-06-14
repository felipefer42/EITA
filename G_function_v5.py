
from time import time
from Interpolate_isochrone import *
import sys
import numpy as np
import scipy.interpolate as intpl
from matplotlib import pyplot as plt
import decimal
import math

__author__ = 'Felipe'

# Version 4.1
# new functions
# load_multiple_isochrones
# interp_multiple_isochrones
# G_for_loaded_isocs

# Version 5.0
# Everything made from zero
# Implementation of isochrones as classes
# Allows the user to freely chose observables

_allowed_colnames_ = ('Z', 'log(age)', 'm_ini', 'logTeff', 'logg', 'Mbol', 'U', 'B', 'V', 'R', 'I')
_allowed_observables_ = ('Z', 'logTeff', 'Teff', 'logg', 'Mbol', 'U', 'B', 'V', 'R', 'I')
_allowed_errors_ = ('Z_err', 'logTeff_err', 'Teff_err', 'logg_err', 'Mbol_err')

Z_sun = 0.0152
Mbol_sun = 4.77
L_sun = 3.846e26  # J/s


def get_filename(t, name_base='set1_res0.01', name_type=1):
    """
    :param t: isochrone age given in years. ex. t=1e9
    :param name_base: base for the name of all isochrones in the set. ex. for isochrone file named 'set3_t5.2e9.dat',
                      name_base should be 'set3'.
    :param name_type: 1 for files names in which age is followed by 'e9' and 2 when 'e9' is omitted.
    :return: file name of isochrone of age t.
    """

    if name_type == 1:
        t_str = "_%0.2f"%(t/1e9) + 'e9'
    elif name_type == 2:
        t_str = "_t%0.1f"%(t/1e9)
    return(name_base + t_str + '.dat')

# Return designated name for an isochrone within an isochrone set
def get_isochrone_name(age):
    """internal use"""
    return str(decimal.Decimal('{:.3f}'.format(age/1e9)).normalize())+'e9'

# Used to check if m_ini is given as one of the columns of an isochrone file
def check_for_m_ini_init(dictionary):
    """Internal use"""
    keys = dictionary.keys()
    if "m_ini" not in keys:
        raise Warning("Column for m_ini not provided must be")


# Used to find the size of each array which is an element of some_array
def list_array_sizes(some_array):
    """internal use"""
    lens = np.zeros(len(some_array))
    for i in range(len(some_array)): lens[i] = len(some_array[i])
    return lens


def Dm_jkl(m):
    """
    Returns array's intervals.
    :param m: array of different masses. Must be ordered.
    :return: difference between successive masses.
    """
    Dm = np.zeros(len(m))
    Dm[1:-1] = m[2:] - m[:-2]
    Dm[0] = 2*(m[1]-m[0])
    Dm[-1] = 2*(m[-1]-m[-2])
    return Dm


def xi(m, a=2.7):
    """
    Kroupa IMF
    :param m: array of masses.
    :param a: alpha parameter for kroupa imf.
    :return: relative distribution of masses.
    """
    return m**(-a)


def abundanceY(Z):
    """
    Return Y abundance for a given Z abundance.
    :param Z: metallicity
    :return: Y abundance scaled for a given metallicity Z.
    """
    return 0.2485 + 1.78*Z


Y_sun = abundanceY(Z_sun)
X_sun = 1-Y_sun-Z_sun


# Considers relation between Y and Z abundances: Y = 0.2485 + 1.78*Z
def MH2Z(MH):
    """
    Transforms metallicity [M/H] to Z
    :param MH: metallicity [M/H] = log(Z/Z_sun) - log(X/X_sun)
    :return: Abundance Z
    """
    k = (Z_sun/X_sun)*10**MH
    Z = (0.7515*k) / (1+2.78*k)
    return Z

# Considers relation between Y and Z abundances: Y = 0.2485 + 1.78*Z
def Z2MH(Z):
    """
    Transforms metallicity Z to [M/H]
    :param Z: metallicity.
    :return: metallicity [M/H] = log(Z/Z_sun) - log(X/X_sun)
    """
    return np.log10( ((Z)/(0.7515-2.78*Z)) * (X_sun/Z_sun) )


# Considers Mbol_sun = 4.77
def fluxbol2magbol(fluxbol, d):
    """
    Transforms bolometric flux to bolometric magnitude
    :param fluxbol: bolometric flux given in units of J/s/m2
    :param d: stellar distance given in m
    :return: stellar bolometric magnitude
    """
    # fluxbol must be given in J/s/m2
    # d must be given in m
    Mbol_sun = 4.77
    return -2.5*np.log10((fluxbol*4*np.pi*d^2)/L_sun) + Mbol_sun


# Considers Mbol_sun = 4.77
def L2magbol(L, unit='L_sun'):
    possible_units = ['L_sun', 'J/s']

    if unit not in possible_units: raise ValueError("unit must be one of the following: {}".format(possible_units))

    Mbol_sun = 4.77

    if unit == 'L_sun':
        return -2.5 * np.log10(L) + Mbol_sun

    if unit == 'J/s':
        L_sun = 3.846e26  # J/s
        return -2.5 * np.log10(L/L_sun) + Mbol_sun


def sample_from_kroupa(N = 10, m0 = 0.01, mf = 10):
    """
    Returns N mass values sampled from a kroupa IMF
    """
    m = np.array([])

    xmax = np.log10(mf)
    xmin = np.log10(m0)

    while len(m) < N:
        x_sort = np.random.random(N)*(xmax-xmin) + xmin  # x = log(m)

        filter_x0 = x_sort < -1.09691  # m < 0.08 M_sun
        filter_x1 = (x_sort >= -1.09691) & (x_sort < -0.30103)  # 0.08 <= m < 0.5 M_Sun
        filter_x2 = x_sort >= -0.30103  # m >= 0.5 M_sun

        y_sort = -np.random.random(N)*4.3  # y = log(Xi(m))

        # Constants chosen to give Xi(m = 0.01) = 1, and maintain continuity
        filter_y0 = y_sort[filter_x0] <= -0.60000 - 0.3 * x_sort[filter_x0]
        filter_y1 = y_sort[filter_x1] <= -1.69691 - 1.3 * x_sort[filter_x1]
        filter_y2 = y_sort[filter_x2] <= -1.99794 - 2.3 * x_sort[filter_x2]

        m = np.concatenate((m, x_sort[filter_x0][filter_y0]))
        m = np.concatenate((m, x_sort[filter_x1][filter_y1]))
        m = np.concatenate((m, x_sort[filter_x2][filter_y2]))

    np.random.shuffle(m)
    m = 10**m

    return m[0:N]

if 0:
    m = sample_from_kroupa(10000)
    print len(m)

def get_nearest(X1, X2):
    """
    Given two arrays X1 and X2, returns an array in which each element corresponds to the element of X2 that is closest
    to the correspondent X1 element.
    """

    Xnew = X1.reshape((len(X1), 1))
    dif = abs(Xnew-X2)

    for i in range(len(X1)):
        Xnew[i] = X2[np.argwhere(dif[i,:] == np.min(dif[i,:]))]

    return Xnew.reshape(len(X1))


def get_argnearest(X1, X2):
    """
    Given two arrays X1 and X2, returns an array in which each element corresponds to the argument of X2 that returns
    the closest value to the correspondent X1 element.
    """

    Xnew = X1.reshape((len(X1), 1))
    dif = abs(Xnew-X2)

    arg = np.zeros(len(Xnew))

    for i in range(len(X1)):
        arg[i] = np.argwhere(dif[i,:] == np.min(dif[i,:]))

    return arg


def load_data(file,
              columns=None,
              type = 'observables',
              **kwargs):
    """
    Returns a dictionary in which keys are the given :param columns.keys() and the values are arrays of observations for
    each key.
    :param     file: path to file containing the that
    :param  columns: dictionary in which keys are observables and values are the correspondent observable column number
                     in :param file.
    :param     type: 'observables' to load observables and 'errors' to load errors.
    :param **kwargs: parameters passed to np.loadtxt.
    """
    colnumbers = columns.values(); colnumbers.sort();
    colnames = []
    for col in colnumbers: colnames.append(columns.keys()[columns.values().index(col)])


    colnumbers = np.array(colnumbers)
    colnames = np.array(colnames)

    if type == 'observables':
        if any(colname not in _allowed_observables_ for colname in colnames):
            raise NameError('Only names in {0} are accepted as colnames'.format(_allowed_observables_))
    elif type == 'errors':
        if any(colname not in _allowed_errors_ for colname in colnames):
            raise NameError('Only names in {0} are accepted as colnames'.format(_allowed_errors_))

    loaded_data = np.loadtxt(file, usecols=colnumbers, **kwargs)
    data_dictionary = {}

    for i in range(len(colnames)):
        data_dictionary[colnames[i]] = loaded_data[:,i]

    return data_dictionary


class Isochrone(object):
    """
    An object containing the information from an isochrone grid for a given age and possible multiple metallicities.
    """

    # TODO implement methods for plotting
    # TODO add more allowed observables
    # TODO method for checking if given age is in agreement with age contained in the file
    # TODO before loading, check if isochrone exists, otherwise, ask if user wants to download it

    def __init__(self,
                 age,
                 path,
                 name_base,
                 name_type=1,
                 columns=None):

        """
        :param       age: isochrone age
        :param      path: path to folder containing isochrone files
        :param name_base: base for the name of all isochrones in the set. ex. for isochrone file named
                          'set3_t5.2e9.dat', name_base should be 'set3'.
        :param name_type: 1 for files names in which age is followed by 'e9' and 2 when 'e9' is omitted.
        :param   columns: dictionary in which keys represent column names and values are the corresponding column number
                          in the isochrone file.

        :return:         self.age = :param age
                       self.Zlist = list of Z values in isochrone's file column 'Z'
                        self.file = {'path': :param path, 'name_base': :param name_base,
                                     'name_type': :param name_type, 'columns': :param columns}
                 self.observables = list of observables listed in columns. ex. 'logTeff', 'logg'
                 for each observable in allowed column names provided in columns:
                    create self.observable

        OBS: allowed column names for the keys of param columns are: 'Z', 'log(age)', 'm_ini', 'logTeff', 'logg', 'Mbol'
        """

        # Dealing with mutable dafaults
        if columns is None:
            columns = {'Z': 0, 'log(age)': 1, 'm_ini': 2, 'logTeff': 5, 'logg': 6, 'Mbol': 7}
        # Warnings
        check_for_m_ini_init(columns)

        # Retrieve colnumbers and colnames
        colnumbers = columns.values(); colnumbers.sort();
        colnames = []
        for col in colnumbers: colnames.append(columns.keys()[columns.values().index(col)])

        colnumbers = np.array(colnumbers)
        colnames = np.array(colnames)

        if any(colname not in _allowed_colnames_ for colname in colnames):
            raise NameError('Only names in {0} are accepted as colnames'.format(_allowed_colnames_))

        Zcol_index = colnumbers[colnames == 'Z'][0]
        Zcol_index = list(colnumbers).index(Zcol_index)

        # Load data
        filename = get_filename(t=age, name_base=name_base, name_type=name_type)
        isoc_data = np.loadtxt(path+'/'+filename, usecols=colnumbers)

        Zcol = isoc_data[:, Zcol_index]
        Zlist = list(set(Zcol)); Zlist.sort()

        # Assign attributes
        self.age = age
        self.Zlist = Zlist
        self.observables = []
        self.file = {'path': path, 'name_base': name_base, 'name_type': name_type, 'columns': columns}

        for i in range(len(colnames)):
            if colnames[i] in _allowed_observables_: self.observables.append(colnames[i])

            attr_temp = []

            for z in Zlist:
                attr_temp.append(isoc_data[Zcol == z,i])

            setattr(self, colnames[i], attr_temp)

        if len(self.observables) < 3:
            raise ValueError('At least three columns of observables must be provided.')


    def __repr__(self):
        Zi = self.Zlist[0]
        Zf = self.Zlist[-1]
        dZ = (self.Zlist[-1]-self.Zlist[0])/(len(self.Zlist)-1)

        return "< Isochrones | age: {:.2f}e9 | Z: from {} to {} by {} | data: {} >".format(self.age/1e9,
                                                                                          Zi, Zf, dZ,
                                                                                          ["m_ini"]+self.observables)


    def logTeff2Teff(self):
        '''
        Transforms isochrone logTeff to Teff
        :return: adds attribute self.Teff
        '''

        if hasattr(self, 'logTeff'):
            self.Teff = []
            for i in range(len(self.Zlist)):
                self.Teff.append(10**self.logTeff[i])
            self.observables.append('Teff')
        else:
            raise AttributeError("Isochrone has no attribute 'logTeff'")


    def plot(self, x_axis='Teff', y_axis='Mbol', fmt='ro', Z=0.012, xlim=None, ylim=None, add=False, invert_x=True,
             invert_y=True, print_Z_msg=True, **kwargs):
        if x_axis not in _allowed_observables_ or y_axis not in _allowed_observables_:
            raise ValueError('parameters x_axis and y_axis must be one of these: {0}'.format(_allowed_observables_))

        #Isochrone Z will be the one closest to the provided :param Z
        Z_isoc = min(self.Zlist, key=lambda x:abs(x-Z))
        Z_isoc_index = self.Zlist.index(Z_isoc)

        if print_Z_msg:
            if Z != Z_isoc:
                print "Isochrone does not contain grid for Z = {0}. Plotting Z = {1} instead.".format(Z, Z_isoc)

        x = getattr(self, x_axis)[Z_isoc_index]
        y = getattr(self, y_axis)[Z_isoc_index]

        plt.plot(x,y,fmt, **kwargs)

        if not add:
            if xlim is None:
                Delta_x = x.max()-x.min()
                xlim = [x.min()-0.04*Delta_x, x.max()+0.04*Delta_x]
            if ylim is None:
                Delta_y = y.max()-y.min()
                ylim = [y.min()-0.04*Delta_y, y.max()+0.04*Delta_y]

            plt.xlim(xlim)
            plt.ylim(ylim)
            if invert_x: plt.gca().invert_xaxis()
            if invert_y: plt.gca().invert_yaxis()
            plt.show()

    def interpolate(self, Ninterp):
        if any(size >= Ninterp for size in list_array_sizes(self.m_ini)): raise ValueError("Ninterp is too small")

        m_ini_interp = []
        for i in range(len(self.m_ini)):
            m_ini_interp.append(np.linspace(self.m_ini[i].min(), self.m_ini[i].max(), Ninterp))

            #Avoid round-offs
            m_ini_interp[i][0] = m_ini_interp[i][0] + (m_ini_interp[i][1]-m_ini_interp[i][0])/10
            m_ini_interp[i][-1] = m_ini_interp[i][-1] - (m_ini_interp[i][-1]-m_ini_interp[i][-2])/10

        for observable in self.observables:
            observable_isoc = getattr(self, observable)
            observable_interp = []
            for i in range(len(observable_isoc)):
                observable_function = intpl.interp1d(self.m_ini[i], observable_isoc[i])

                observable_interp.append(observable_function(m_ini_interp[i]))

            setattr(self, observable, observable_interp)
        setattr(self, 'm_ini', m_ini_interp)

    def reset_interpolation(self):
        self.__init__(age=self.age,
                      path=self.file['path'],
                      name_base=self.file['name_base'],
                      name_type = self.file['name_type'],
                      columns=self.file['columns'])

        if hasattr(self, 'Teff'):
            self.logTeff2Teff()


    def Gj_interp(self, observations, errors, Ninterp = 10000):
        # Observations must be a dictionary
        if any(observation not in _allowed_observables_ for observation in observations.keys()):
            raise NameError('Only names in {0} are accepted as observations'.format(_allowed_observables_))

        if any(error not in _allowed_errors_ for error in errors.keys()):
            raise NameError('Only names in {0} are accepted as errors'.format(_allowed_errors_))

        # Perform interpolation
        m_ini_interp = []
        for i in range(len(self.m_ini)):
            m_ini_interp.append(np.linspace(self.m_ini[i].min(), self.m_ini[i].max(), Ninterp))
            #Avoid round-offs
            m_ini_interp[i][0] = m_ini_interp[i][0] + (m_ini_interp[i][1]-m_ini_interp[i][0])/10
            m_ini_interp[i][-1] = m_ini_interp[i][-1] - (m_ini_interp[i][-1]-m_ini_interp[i][-2])/10

        observables_interp = {}

        for observable in self.observables:
            observable_isoc = getattr(self, observable)
            observable_interp = []
            for i in range(len(observable_isoc)):
                observable_function = intpl.interp1d(self.m_ini[i], observable_isoc[i])
                observable_interp.append(observable_function(m_ini_interp[i]))
            observables_interp[observable] = observable_interp

        # Create dm arrays
        dm = []
        for i in range(len(m_ini_interp)):
            dm.append(Dm_jkl(m_ini_interp[i]))

        # Concatenate different metallicities
        dm_isoc = np.array([])
        m_ini_isoc = np.array([])

        observables_isoc = {}

        for observable in observations.keys():
            observables_isoc[observable] = np.array([])

        for i in range(len(m_ini_interp)):
            dm_isoc = np.concatenate((dm_isoc, dm[i]))
            m_ini_isoc = np.concatenate((m_ini_isoc, m_ini_interp[i]))
            for observable in observations.keys():
                observables_isoc[observable] = np.concatenate((observables_isoc[observable],
                                                               observables_interp[observable][i]))

        if type (observations[observations.keys()[0]]) == int or type(observations[observations.keys()[0]]) == float:
            N_lines = 1
        else:
            N_lines = len(observations[observations.keys()[0]])
        N_cols = len(m_ini_isoc)

        for key in observations.keys(): observations[key] = np.array(observations[key])  # for single object problem
        for key in errors.keys(): errors[key] = np.array(errors[key])  # for single object problem

        # Reshape observations
        for observable in observations.keys():
            observations[observable] = observations[observable].reshape((N_lines,1))
            errors[observable+'_err'] = errors[observable+'_err'].reshape((N_lines,1))
        Chi2 = np.zeros((N_lines, N_cols))

        # Calculate Xsi2
        for observable in observations.keys():
            Chi_arg = (observables_isoc[observable] - observations[observable]) / errors[observable+'_err']
            Chi_arg = Chi_arg*Chi_arg
            Chi2 = Chi2 + Chi_arg

        #Calculate L
        L = np.exp(-Chi2/2)
        sqrt2pi = np.sqrt(2*np.pi)
        for observable in observations.keys():
            L = (1 / (sqrt2pi*errors[observable+'_err'])) * L

        #Calculate Gj
        Gj = L*xi(m_ini_isoc)*dm_isoc
        Gj = Gj.sum(1).reshape((N_lines, 1))

        return Gj


    def Gj(self, observations, errors):
        # Observations must be a dictionary
        if any(observation not in _allowed_observables_ for observation in observations.keys()):
            raise NameError('Only names in {0} are accepted as observations'.format(_allowed_observables_))

        if any(error not in _allowed_errors_ for error in errors.keys()):
            raise NameError('Only names in {0} are accepted as errors'.format(_allowed_errors_))

        # Create dm arrays
        dm = []
        for i in range(len(self.m_ini)):
            dm.append(Dm_jkl(self.m_ini[i]))

        # Concatenate different metallicities
        dm_isoc = np.array([])
        m_ini_isoc = np.array([])

        observables_isoc = {}

        for observable in observations.keys():
            observables_isoc[observable] = np.array([])

        for i in range(len(self.m_ini)):
            dm_isoc = np.concatenate((dm_isoc, dm[i]))
            m_ini_isoc = np.concatenate((m_ini_isoc, self.m_ini[i]))
            for observable in observations.keys():
                observable_i = getattr(self, observable)[i]
                observables_isoc[observable] = np.concatenate((observables_isoc[observable], observable_i))


        if type (observations[observations.keys()[0]]) == int or type(observations[observations.keys()[0]]) == float:
            N_lines = 1
        else:
            N_lines = len(observations[observations.keys()[0]])
        N_cols = len(m_ini_isoc)

        for key in observations.keys(): observations[key] = np.array(observations[key])  # for single object problem
        for key in errors.keys(): errors[key] = np.array(errors[key])  # for single object problem

        # Reshape observations
        for observable in observations.keys():
            observations[observable] = observations[observable].reshape((N_lines,1))
            errors[observable+'_err'] = errors[observable+'_err'].reshape((N_lines,1))
        Chi2 = np.zeros((N_lines, N_cols))

        # Calculate Xsi2
        for observable in observations.keys():
            Chi_arg = (observables_isoc[observable] - observations[observable]) / errors[observable+'_err']
            Chi_arg = Chi_arg*Chi_arg
            Chi2 = Chi2 + Chi_arg

        #Calculate L
        L = np.exp(-Chi2/2)
        sqrt2pi = np.sqrt(2*np.pi)
        for observable in observations.keys():
            L = (1 / (sqrt2pi*errors[observable+'_err'])) * L

        #Calculate Gj
        Gj = L*xi(m_ini_isoc)*dm_isoc
        Gj = Gj.sum(1).reshape((N_lines, 1))

        return Gj


class Isochrone_Set(object):

    # TODO implement methods for plotting

    def __init__(self,
                 age0,
                 agef,
                 dage,
                 path,
                 name_base,
                 name_type,
                 columns = None,
                 display_time = True,
                 print_steps = True):

        if display_time: t0 = time()

        if columns is None:
            columns = {'Z': 0, 'log(age)': 1, 'm_ini': 2, 'logTeff': 5, 'logg': 6, 'Mbol': 7}

        self.isochrones = {}
        self.ages = np.arange(age0, agef, dage)

        if print_steps:
            print "Loading isochrones."
            N_isocs = len(self.ages)
            N_load = 0.

        for age in self.ages:
            if print_steps:
                pct = (N_load/N_isocs)*100
                sys.stdout.write('\r%3d%%'%(pct))
                sys.stdout.flush()
                N_load += 1

            isochrone_name = str(decimal.Decimal('{:.3f}'.format(age/1e9)).normalize())+'e9'
            self.isochrones[isochrone_name] = Isochrone(age=age,
                                                        path=path,
                                                        name_base=name_base,
                                                        name_type=name_type,
                                                        columns=columns)

        if print_steps:
            pct = 100
            sys.stdout.write('\r%3d%%\n'%(pct))
            sys.stdout.flush()

        if display_time: tf = time(); print('Loading process took {0} sec.\n'.format(tf-t0))


    def interpolate(self,
                    Ninterp = 10000,
                    display_time = True,
                    print_steps = True):

        if display_time: t0 = time()

        if print_steps:
            print "Interpolating isochrones."
            N_isocs = len(self.ages)
            N_load = 0.

        for age in self.ages:
            if print_steps:
                pct = (N_load/N_isocs)*100
                sys.stdout.write('\r%3d%%'%(pct))
                sys.stdout.flush()
                N_load += 1

            isochrone_name = str(decimal.Decimal('{:.3f}'.format(age/1e9)).normalize())+'e9'
            try: self.isochrones[isochrone_name].interpolate(Ninterp)
            except ValueError:
                for age_reset in self.ages:
                    isochrone_name_reset = str(decimal.Decimal('{:.3f}'.format(age_reset/1e9)).normalize())+'e9'
                    self.isochrones[isochrone_name_reset].reset_interpolation()
                raise ValueError("Ninterp may be too small")

        if print_steps:
            pct = 100
            sys.stdout.write('\r%3d%%\n'%(pct))
            sys.stdout.flush()

        if display_time: tf = time(); print("Interpolation process took {0} sec.\n".format(tf-t0))

    def __repr__(self):
        first_isochrone_name = str(decimal.Decimal('{:.3f}'.format(self.ages[0]/1e9)).normalize())+'e9'
        first_isochrone = self.isochrones[first_isochrone_name]
        Zi = first_isochrone.Zlist[0]
        Zf = first_isochrone.Zlist[-1]
        dZ = (first_isochrone.Zlist[-1]-first_isochrone.Zlist[0])/(len(first_isochrone.Zlist)-1)
        observables = first_isochrone.observables

        agei = self.ages[0]/1e9
        agef = self.ages[-1]/1e9
        dage = ((self.ages[-1]-self.ages[0])/(len(self.ages)-1))/1e9
        return "< Isochrones | age (yr): from {:.2f}e9 to {:.2f}e9 by {:.2f}e9 |" \
               " Z: from {:.4f} to {:.4f} by {:.4f} | data: {} >".format(agei, agef, dage,
                                                                         Zi, Zf, dZ,
                                                                         ["m_ini"]+observables)

    def plot(self, x_axis='Teff', y_axis='Mbol', fmt='b-', Z=0.012, xlim=None, ylim=None, add=False, invert_x=True,
             invert_y=True, print_Z_msg=True, **kwargs):

        # First isochrone:
        age = self.ages[0]
        isochrone = self.isochrones[get_isochrone_name(age)]

        if xlim is None:
            # Get x range for this age
            xmin_list = []
            xmax_list = []
            for xlist in getattr(isochrone, x_axis):
                xmin_list.append(xlist.min())
                xmax_list.append(xlist.max())

            xmin = np.array(xmin_list).min()
            xmax = np.array(xmax_list).max()

        if ylim is None:
            # Get y range for this age
            ymin_list = []
            ymax_list = []
            for ylist in getattr(isochrone, y_axis):
                ymin_list.append(ylist.min())
                ymax_list.append(ylist.max())

            ymin = np.array(ymin_list).min()
            ymax = np.array(ymax_list).max()

        # Plot first isochrone
        isochrone.plot(x_axis=x_axis, y_axis=y_axis, fmt=fmt, Z=Z, add=True, print_Z_msg=print_Z_msg, **kwargs)

        for i in range(1, len(self.ages)):
            age = self.ages[i]
            isochrone = self.isochrones[get_isochrone_name(age)]

            if xlim is None:
                # Get x range for this age
                xmin_list = []
                xmax_list = []
                for xlist in getattr(isochrone, x_axis):
                    xmin_list.append(xlist.min())
                    xmax_list.append(xlist.max())

                xmin_list = np.array(xmin_list)
                xmax_list = np.array(xmax_list)

                # Update x range for all ages
                if xmin_list.min() < xmin: xmin = xmin_list.min()
                if xmax_list.max() > xmax: xmax = xmax_list.max()

            if ylim is None:
                # Get y range for this age
                ymin_list = []
                ymax_list = []
                for ylist in getattr(isochrone, y_axis):
                    ymin_list.append(ylist.min())
                    ymax_list.append(ylist.max())

                ymin_list = np.array(ymin_list)
                ymax_list = np.array(ymax_list)

                # Update x range for all ages
                if ymin_list.min() < ymin: ymin = ymin_list.min()
                if ymax_list.max() > ymax: ymax = ymax_list.max()

            # Plot isochrone for this age
            isochrone.plot(x_axis=x_axis, y_axis=y_axis, fmt=fmt, Z=Z, add=True, print_Z_msg=False, **kwargs)

        if not add:
            # Finish plot
            plt.xlim(xlim)
            plt.ylim(ylim)
            if invert_x: plt.gca().invert_xaxis()
            if invert_y: plt.gca().invert_yaxis()
            plt.show()

    def logTeff2Teff(self):
        for age in self.ages:
            self.isochrones[get_isochrone_name(age)].logTeff2Teff()

    def sample(self,  N=1000, agedist = 'uniform', Zdist = 'uniform', imf = 'kroupa', m_min = -np.inf, m_max = np.inf,
               age_min = 0.1e9, age_max = 13.5e9, Z_min = 0.0001, Z_max = 0.0401, print_steps = True,
               display_time = True):

        if display_time: t0 = time()

        if print_steps:
            sys.stdout.write('Sampling {0} stars from isochrone set\n'.format(N))
            sys.stdout.flush()

        Zlist = np.array(self.isochrones[self.isochrones.keys()[0]].Zlist)

        age0 = age_min if age_min >= self.ages.min() else self.ages.min()
        agef = age_max if age_max <= self.ages.max() else self.ages.max()

        Z0 = Z_min if Z_min >= Zlist.min() else Zlist.min()
        Zf = Z_max if Z_max <= Zlist.max() else Zlist.max()

        if agedist == 'uniform':
            t_sampled = np.random.random(N)*(agef-age0) + age0

            if Zdist == 'uniform':
                Z_sampled = np.random.random(N)*(Zf-Z0) + Z0

        t_sampled = get_nearest(t_sampled, self.ages)
        Z_sampled = get_nearest(Z_sampled, Zlist)

        sampled = {} # Dictionary that will receive sampled data
        sampled['age'] = t_sampled
        sampled['Z'] = Z_sampled
        sampled['m_ini'] = np.array([])

        isoc_name = get_isochrone_name(self.ages[0])
        isochrone = self.isochrones[isoc_name]
        for observable in isochrone.observables:
            if observable != 'Z':
                sampled[observable] = np.array([])

        for i in range(N):
            isoc_name = get_isochrone_name(sampled['age'][i])
            isochrone = self.isochrones[isoc_name]

            Z = sampled['Z'][i]
            arg_Z = int(np.arange(len(Zlist))[Zlist == Z])
            if Zlist[arg_Z] != Z: print "deu treta ai"

            m_isoc = isochrone.m_ini[arg_Z]
            m0 = m_isoc.min() if m_isoc.min() > m_min else m_min
            mf = m_isoc.max() if m_isoc.max() < m_max else m_max

            while m0 >= m_isoc.max(): m0 -= 5*m_isoc.max()/100.0
            while mf <= m_isoc.min(): mf += 5*m_isoc.min()/100.0

            m_sampled_i = sample_from_kroupa(1, m0, mf)
            sampled['m_ini'] = np.concatenate((sampled['m_ini'], m_sampled_i))

            for observable in isochrone.observables:
                if observable != 'Z':
                    observable_isoc = getattr(isochrone, observable)[arg_Z]
                    observable_function = intpl.interp1d(m_isoc, observable_isoc)
                    sampled[observable] = np.concatenate((sampled[observable], observable_function(m_sampled_i)))

            if print_steps:
                sys.stdout.write("\rSampled {0} stars of {1}.".format(i, N))
                sys.stdout.flush()

        if print_steps:
            sys.stdout.write("\rSampled {0} stars of {1}.".format(i+1, N))
            sys.stdout.flush()

        if display_time: tf = time(); print("\nSampling {0} stars took {1} sec.\n".format(N, tf-t0))

        return sampled


    def G(self,
          observations,
          errors,
          display_time = True,
          print_steps = True,
          interp = False,
          Ninterp = 10000,
          N_obj_max = 1000):

        if display_time: t0 = time()


        N_isocs = len(self.ages)

        if type (observations[observations.keys()[0]]) == int or type(observations[observations.keys()[0]]) == float:
            N_obj = 1
        else:
            N_obj = len(observations[observations.keys()[0]])

        for key in observations.keys(): observations[key] = np.array(observations[key])  # for single object problem
        for key in errors.keys(): errors[key] = np.array(errors[key])  # for single object problem

        if print_steps: print "Calculating G(t) for {0} stars.".format(N_obj)

        G = np.zeros((N_obj, N_isocs))

        if N_obj_max < N_obj:
            N_steps = int(math.ceil(float(N_obj)/N_obj_max))
        else:
            N_steps = 1

        for step in range(N_steps):
            obj_0 = step*N_obj_max
            if step == (N_steps-1):
                obj_f = N_obj
            else:
                obj_f = (step+1)*N_obj_max

            N_obj_step = obj_f - obj_0

            if print_steps: print "Calculating G(t) for objects {0} to {1}.".format(obj_0, obj_f)

            observations_step = dict(observations)
            errors_step = dict(errors)

            for observation in observations_step.keys():
                observations_step[observation] = observations_step[observation][obj_0:obj_f]
                errors_step[observation+'_err'] = errors_step[observation+'_err'][obj_0:obj_f]

            for i in range(len(self.ages)):
                if print_steps:
                    pct = (float(i)/N_isocs)*100
                    sys.stdout.write('\r%3d%%'%(pct))
                    sys.stdout.flush()

                age = self.ages[i]
                isochrone_name = str(decimal.Decimal('{:.3f}'.format(age/1e9)).normalize())+'e9'

                if interp == True:
                    Gj = self.isochrones[isochrone_name].Gj_interp(observations=observations_step,
                                                                   errors=errors_step,
                                                                   Ninterp=Ninterp)
                else:
                    Gj = self.isochrones[isochrone_name].Gj(observations=observations_step,
                                                            errors=errors_step)

                Gj = Gj.reshape((1, N_obj_step))
                G[obj_0:obj_f,i] = Gj

            if print_steps:
                pct = 100
                sys.stdout.write('\r%3d%%\n'%(pct))
                sys.stdout.flush()

        if print_steps: tf = time(); print "Calculating G(t) for {0} stars took {1} sec.\n".format(N_obj, tf-t0)

        return G


######################################################################################################

# Test isochrone.sample:
if 0:
    isocs3 = Isochrone_Set(age0=0.1e9, agef=13.6e9, dage=0.1e9,
                      path='/home/felipe/Documents/Isochrones/Sets/Set3',
                      name_base='Parsec_Set3', name_type=2)

    isocs.logTeff2Teff()

    sample = isocs.sample(N = 50, Z_min=0.02, Z_max = 0.022, m_min = 0.9, m_max = 1.5)
    plt.plot(sample['Teff'], sample['Mbol'], 'ro')
    isocs.plot(Z = 0.02)

    isocs = Isochrone_Set(age0=0.1e9, agef=13.6e9, dage=0.1e9,
                          path='/home/felipe/Documents/Isochrones/Sets/Set3',
                          name_base='Parsec_Set3', name_type=2)
    isocs.logTeff2Teff()
    observations = {'Teff': sample['Teff'], 'Mbol': sample['Mbol'], 'Z': sample['Z']}
    Teff_err = [50]*50
    Mbol_err = [0.1]*50
    Z_err = [0.001]*50
    errors = {'Teff_err': Teff_err, 'Mbol_err': Mbol_err, 'Z_err': Z_err}

    isocs.interpolate()
    G = isocs.G(observations, errors, interp = False)

    tseq = np.arange(0.1e9, 13.6e9, 0.1e9)
    tML = []
    tE = []
    for i in range(10):
        tML.append(tseq[G[i,:] == max(G[i,:])][0])
        norm = sum(G[i,:])
        tE.append(sum(tseq*G[i,:])/norm)

    i = 9

    print sample['age'][i]/1e9
    print tML[i]/1e9
    print tE[i]/1e9
    plt.plot(tseq, G[i,:], 'g-')
    plt.show()


# Test without mbol
if 0:
    observations = load_data(file="/home/felipe/Documents/Isochrones/Teste G_function_v5/AmostraApronta20.tsv",
                             columns={'Teff': 1, 'logg': 3, 'Z': 5, 'Mbol': 7},
                             type = 'observables',
                             delimiter='|')

    errors = load_data(file="/home/felipe/Documents/Isochrones/Teste G_function_v5/AmostraApronta20.tsv",
                       columns={'Teff_err': 2, 'logg_err': 4, 'Z_err': 6, 'Mbol_err': 8},
                       type = 'errors',
                       delimiter="|")

    # Plot HR diagram
    #isocs = Isochrone_Set(age0=0.1e9, agef=13.5e9, dage=1e9,
    #                      path='/home/felipe/Documents/Isochrones/Sets/Set3',
    #                      name_base='Parsec_Set3', name_type=2)


    #isocs.logTeff2Teff()
    #plt.plot(observations['Teff'], observations['logg'], 'ro')
    #isocs.plot('Teff', 'logg', fmt = 'b-', Z = 0.0152,  xlim=[4500,8500], ylim=[0,7])

    # Load Isochrones
    isocs = Isochrone_Set(age0=0.1e9, agef=13.5e9, dage=0.1e9,
                          path='/home/felipe/Documents/Isochrones/Sets/Set3',
                          name_base='Parsec_Set3', name_type=2)
    isocs.logTeff2Teff()

    G_interp = isocs.G(observations=observations, errors=errors, interp=True, N_obj_max = 10)

    # Interpolate Isochrones
    isocs.interpolate(Ninterp=10000)
    # Calculate G
    G = isocs.G(observations=observations, errors=errors, N_obj_max = 10)
    G2 = isocs.G(observations=observations, errors=errors, N_obj_max = 20)

    print G2 - G_interp
    print G
    # Save Results
    #np.save('/home/felipe/Documents/Isochrones/Teste G_function_v5/G_GCS_200_no_mbol.npy', G)
    #np.save('/home/felipe/Documents/Isochrones/Teste G_function_v5/G_GCS_200_no_mbol.npy', G)

# Test create Isochrone object
if 0:
    isoc = Isochrone(age=2e9,
                     path='/home/felipe/Documents/Isochrones/Sets/Set3',
                     name_base='Parsec_Set3',
                     name_type=2,
                     columns={'Z': 0, 'logTeff': 5, 'logg': 6, 'm_ini':2, 'Mbol':7})

    print isoc.observables
    print len(isoc.logTeff[0])
    #isoc.interpolate(1000)
    print len(isoc.logTeff[0])

    observations = {'Z': 0.012, 'logTeff': 3.5, 'Mbol': 4}
    errors = {'Z_err': 0.001, 'logTeff_err': 0.1, 'Mbol_err':0.1}

    Gj = isoc.Gj(observations, errors)
    print Gj

# Test plot single isochrone
if 0:
    isoc2 = Isochrone(age=2e9,
                     path='/home/felipe/Documents/Isochrones/Sets/Set3',
                     name_base='Parsec_Set3',
                     name_type=2,
                     columns={'Z': 0, 'logTeff': 5, 'logg': 6, 'm_ini':2, 'Mbol':7})

    isoc2.logTeff2Teff()
    isoc2.interpolate(10000)
    isoc2.plot('Teff', 'Mbol', 'r-', add = True)

    isoc = Isochrone(age=1e9,
                     path='/home/felipe/Documents/Isochrones/Sets/Set3',
                     name_base='Parsec_Set3',
                     name_type=2,
                     columns={'Z': 0, 'logTeff': 5, 'logg': 6, 'm_ini':2, 'Mbol':7})


    isoc.logTeff2Teff()
    isoc.interpolate(10000)
    isoc.plot('Teff', 'Mbol', 'b-')

# Test plot of isochrone set
if 0:
    isocs = Isochrone_Set(age0=0.1e9, agef=13.5e9, dage=0.1e9,
                      path='/home/felipe/Documents/Isochrones/Sets/Set3',
                      name_base='Parsec_Set3', name_type=2)

    print isocs
    isocs.logTeff2Teff()
    #isocs.plot('Teff', 'logg', add = True, fmt = 'g-', Z = 0.001)
    #isocs.plot('Teff', 'logg', Z = 0.01)
    sample = isocs.sample(1000)

    plt.plot(sample['Teff'], sample['Mbol'], 'ro')

    isocs_plot = Isochrone_Set(age0=0.5e9, agef=13.6e9, dage=2e9,
                      path='/home/felipe/Documents/Isochrones/Sets/Set3',
                      name_base='Parsec_Set3', name_type=2)
    isocs_plot.logTeff2Teff()
    isocs_plot.interpolate(Ninterp=2000)
    isocs_plot.plot('Teff', 'Mbol', add = True, fmt = 'b-', Z = 0.0001)
    isocs_plot.plot('Teff', 'Mbol', add = True, fmt = 'g-', Z = 0.0150)
    #isocs.plot('Teff', 'Mbol', fmt = 'c-', xlim=[2000,15000], ylim=[0.75,6.75], Z = 0.03)
    isocs_plot.plot('Teff', 'Mbol', fmt = 'c-', Z = 0.03)
    plt.show()


if 0:
    isocs = Isochrone_Set(age0=0.1e9, agef=13.4e9, dage=0.1e9,
                      path='/home/felipe/Documents/Isochrones/Sets/Set3',
                      name_base='Parsec_Set3', name_type=2)

    isocs.interpolate(Ninterp=10000)

    observations = {'Z': np.array([0.012, 0.015, 0.018]), 'logTeff': np.log10(np.array([4000, 4500, 5000])),
                    'Mbol': np.array([4, 4.1, 4.2])}
    errors = {'Z_err': np.array([0.01, 0.01, 0.01]), 'logTeff_err': np.array([0.1, 0.1, 0.1]),
              'Mbol_err': np.array((0.1, 0.1, 0.1))}

    G = isocs.G(observations, errors)
    print G

def load_isochrone(t,
                   path,
                   name_base = 'set1_res0.01',
                   usecols=(0,1,2,5,7),
                   nametype = 1
                   ):



    filename = get_filename(t, name_base, name_type=nametype)
    data = np.loadtxt(path+'/'+filename, usecols=usecols)

    Z = data[:,0]
    m = data[:,2]
    Mv = data[:,4]
    Te = 10**data[:,3]

    # Sorting data according to mass (it's important when evaluating discrete integral that uses dm)
    sort = np.argsort(m)

    Z = Z[sort]
    m = m[sort]
    Mv = Mv[sort]
    Te = Te[sort]

    return(Z, m, Mv, Te)

#######################
# Prior distributions #
#######################
#                     #
def xi(m, a = 2.7):
    return m**(-a)

def phi(Z):
    return 1/Z

def psi(t):
    return 1
#                     #
#######################

########################################################################
# G_j                                                                  #
########################################################################
#                                                                      #
def Xsi2_jkl(Te, Te_obs, Te_err, Mv, Mv_obs, Mv_err, Z, Z_obs, Z_err):
    arg_Te = ((Te_obs-Te)/Te_err)**2
    arg_Mv = ((Mv_obs-Mv)/Mv_err)**2
    arg_Z =((Z_obs-Z)/Z_err)**2

    return(arg_Te+arg_Mv+arg_Z)

def L_jkl(Te, Te_obs, Te_err, Mv, Mv_obs, Mv_err, Z, Z_obs, Z_err):
    sqrt2pi = np.sqrt(2*np.pi)

    arg_Te_err = 1/(sqrt2pi*Te_err)
    arg_Mv_err = 1/(sqrt2pi*Mv_err)
    arg_Z_err = 1/(sqrt2pi*Z_err)

    arg_err = arg_Te_err*arg_Mv_err*arg_Z_err
    Xsi2_jkl_ = Xsi2_jkl(Te=Te, Te_obs=Te_obs, Te_err=Te_err,
                         Mv=Mv, Mv_obs=Mv_obs, Mv_err=Mv_err,
                         Z=Z, Z_obs=Z_obs, Z_err=Z_err)

    return(arg_err*np.exp(-Xsi2_jkl_/2))

def Dm_jkl(m):
    Dm = np.zeros(len(m))
    Dm[1:-1] = m[2:] - m[:-2]
    Dm[0] = 2*(m[1]-m[0])
    Dm[-1] = 2*(m[-1]-m[-2])
    return Dm

def G_j(Te, Te_obs, Te_err, Mv, Mv_obs, Mv_err, Z, Z_obs, Z_err, m, a = 2.7):
    L_jkl_ = L_jkl(Te=Te, Te_obs=Te_obs, Te_err=Te_err,
                 Mv=Mv, Mv_obs=Mv_obs, Mv_err=Mv_err,
                 Z=Z, Z_obs=Z_obs, Z_err=Z_err)

    xi_jkl_ = xi(m, a)
    phi_k_ = phi(Z)

    Dm_jkl_ = Dm_jkl(m)
    G_jkl_ =  L_jkl_ * xi_jkl_ * phi_k_ * Dm_jkl_
    G_j_ = G_jkl_.sum()

    return(G_j_)
#                                                                      #
########################################################################

#################################################################################################
# G(t)                                                                                          #
#################################################################################################
#                                                                                               #
def G(Te_obs, Te_err, Mv_obs, Mv_err, Z_obs, Z_err,
      path, name_base, ti = 0.1e9, tf = 13.5e9, dt = 0.1e9,
      nsigma = 20, nintpl = 10000, a = 2.7,
      plot = True, get_time = True):

    if get_time: t0 = time()
    ############################################################
    # Add later                                                #
    # Range for which isochrones values are considered
    #Z_lim = [Z_obs-nsigma*Z_err, Z_obs + nsigma*Z_err]
    #Te_lim = [Te_obs - nsigma*Te_err, Te_obs + nsigma*Te_err]
    #Mv_lim = [Mv_obs - nsigma*Mv_err, Mv_obs + nsigma*Mv_err]
    #                                                          #
    ############################################################

    # Times t for which G(t) will be calculated
    t_seq = np.arange(ti, tf, dt)

    # Vector that will receive values of G(t)
    G = np.zeros(len(t_seq))

    # For each age tj
    for j in range(len(t_seq)):
        t_j = t_seq[j]
        print('Working with %0.2f'%(t_j/1e9)+'e9 Ga isochrone.')

        # Loading isochrone tj
        Z, m, Mv, Te = load_isochrone(t_j, path, name_base)

        # Isochrones metallicity range
        Z_seq = list(set(Z))

        #####################
        # Add later         #
        # Consider Filters  #
        #####################

        # For each metallicity Zk
        for k in range(len(Z_seq)):
            Zk = Z_seq[k]

            # Filtering isochrone metallicity value
            filter_Zk = Z == Zk
            Te_filt = Te[filter_Zk]
            Mv_filt = Mv[filter_Zk]
            m_filt = m[filter_Zk]

            # Interpolating values
            Te_intpl, Mv_intpl, m_intpl = interp_isoc(Te_filt, Mv_filt, m_filt, nintpl = nintpl)

            # Calculating Gjk
            G[j] += G_j(Te=Te_intpl, Te_obs=Te_obs, Te_err=Te_err,
                        Mv=Mv_intpl, Mv_obs=Mv_obs, Mv_err=Mv_err,
                        Z=Zk, Z_obs=Z_obs, Z_err=Z_err,
                        m=m_intpl, a=a)

        print('Obtained G(j) = {0}'.format(G[j]))

    # Normalizing
    G = G/G.max()

    if get_time: tf = time()

    # Plot result
    if plot:
        plt.plot(t_seq, G)
        plt.show()

    if get_time: print("It took = {0} sec.".format(tf-t0))

    return G
#                                                                                               #
#################################################################################################

#################################################################################################
# Load all isochrones                                                                           #
#################################################################################################
#                                                                                               #
def load_multiple_isochrones(tseq, path, name_base, nametype):

    N_isocs = len(tseq)

    Z_isocs = []
    m_isocs = []
    Mv_isocs = []
    Te_isocs = []

    #Loading isochrones
    for i in range(N_isocs):
        pct = (float(i)/N_isocs)*100
        sys.stdout.write('\r%3d%%'%(pct))
        sys.stdout.flush()

        Z_isoc_i, m_isoc_i, Mv_isoc_i, Te_isoc_i = load_isochrone(tseq[i], path, name_base, nametype=nametype)
        Z_isocs.append(Z_isoc_i)
        m_isocs.append(m_isoc_i)
        Mv_isocs.append(Mv_isoc_i)
        Te_isocs.append(Te_isoc_i)

    return (Z_isocs, m_isocs, Mv_isocs, Te_isocs)
#                                                                                               #
#################################################################################################

#################################################################################################
# Interpolate multiple isochrones                                                               #
#################################################################################################
#                                                                                               #
def interp_multiple_isochrones(Z_isocs, m_isocs, Mv_isocs, Te_isocs, nintpl = 10000):

    N_isocs = len(Z_isocs)

    Z_interp = []
    m_interp = []
    Mv_interp = []
    Te_interp = []

    for i in range(N_isocs):
        pct = (float(i)/N_isocs)*100
        sys.stdout.write('\r%3d%%'%(pct))
        sys.stdout.flush()

        Z_interp_i = []
        m_interp_i = []
        Mv_interp_i = []
        Te_interp_i = []

        Z_isoc_i = Z_isocs[i]
        m_isoc_i = m_isocs[i]
        Mv_isoc_i = Mv_isocs[i]
        Te_isoc_i = Te_isocs[i]

        Zseq_i = list(set(Z_isoc_i))
        for Zk in Zseq_i:
            filter_Zk = Z_isoc_i == Zk

            m_isoc_i_filt = m_isoc_i[filter_Zk]
            Mv_isoc_i_filt = Mv_isoc_i[filter_Zk]
            Te_isoc_i_filt = Te_isoc_i[filter_Zk]

            Mv_func = intpl.interp1d(m_isoc_i_filt, Mv_isoc_i_filt)
            Te_func = intpl.interp1d(m_isoc_i_filt, Te_isoc_i_filt)

            m_intpl = np.linspace(m_isoc_i_filt[0], m_isoc_i_filt[-1], nintpl)
            Mv_interp_i = Mv_interp_i + list(Mv_func(m_intpl))
            Te_interp_i = Te_interp_i + list(Te_func(m_intpl))
            Z_interp_i = Z_interp_i + [Zk]*nintpl
            m_interp_i = m_interp_i + list(m_intpl)

        Z_interp_i = np.array(Z_interp_i)
        Mv_interp_i = np.array(Mv_interp_i)
        Te_interp_i = np.array(Te_interp_i)
        m_interp_i = np.array(m_interp_i)

        Z_interp.append(Z_interp_i)
        Mv_interp.append(Mv_interp_i)
        Te_interp.append(Te_interp_i)
        m_interp.append(m_interp_i)

    return(Z_interp, m_interp, Mv_interp, Te_interp)
#                                                                                               #
#################################################################################################

if 0:
    t0_total = time()
    path=r"C:\Users\Felipe\Documents\Doutorado\Projeto Idade Isocronal\Isochrone\Isochrones_set_parsec1.0_res0.1Ga"
    name_base = "parsec1.0_res0.5"

    tseq = np.arange(0.1e9, 13.5e9, 0.1e9)

    print('Loading isochrones')
    t0_load = time()
    Z_isocs, m_isocs, Mv_isocs, Te_isocs = load_multiple_isochrones(tseq, path, name_base)
    tf_load = time()
    print('Loading process took {0} sec.'.format(tf_load-t0_load))

    print('Interpolating isochrones')
    t0_intpl = time()
    Z_interp, m_interp, Mv_interp, Te_interp = interp_multiple_isochrones(Z_isocs, m_isocs, Mv_isocs, Te_isocs, nintpl = 10000)
    tf_intpl = time()
    print('Interpolation process took {0} sec.'.format(tf_intpl-t0_intpl))

    tf_total = time()
    print('Total process took {0} sec.'.format(tf_total-t0_total))

#################################################################################################
# G(t) with already loaded isochrones                                                           #
#################################################################################################
#
def G_for_loaded_isocs(Te_obs, Te_err, Mv_obs, Mv_err, Z_obs, Z_err,
                       Te_isocs, Mv_isocs, Z_isocs, m_isocs,
                       tseq = np.arange(0.1, 13.5, 0.1), a = 2.7,
                       get_time = True):

    if get_time: t0 = time()
    G = np.zeros(len(tseq))
    for j in range(len(tseq)):
        tj = tseq[j]

        Z = Z_isocs[j]
        Mv = Mv_isocs[j]
        Te = Te_isocs[j]
        m = m_isocs[j]

        G[j] = G_j(Te=Te, Te_obs=Te_obs, Te_err=Te_err,
                   Mv=Mv, Mv_obs=Mv_obs, Mv_err=Mv_err,
                   Z=Z, Z_obs=Z_obs, Z_err=Z_err,m=m,a=a)

    #Normalizing
    G = G/G.max()

    if get_time:
        tf = time()
        print('It took {0} sec. to calculate G(j)'.format(tf-t0))

    return G
#                                                                                               #
#################################################################################################

if 0:
    t0_total = time()
    path=r"C:\Users\Felipe\Documents\Doutorado\Projeto Idade Isocronal\Isochrone\Isochrones_set_parsec1.0_res0.1Ga"
    name_base = "parsec1.0_res0.5"

    tseq = np.arange(0.1e9, 13.5e9, 0.1e9)

    print('Loading isochrones')
    t0_load = time()
    Z_isocs, m_isocs, Mv_isocs, Te_isocs = load_multiple_isochrones(tseq, path, name_base)
    tf_load = time()
    print('Loading process took {0} sec.\n'.format(tf_load-t0_load))

    print('Interpolating isochrones')
    t0_intpl = time()
    Z_isocs, m_isocs, Mv_isocs, Te_isocs = interp_multiple_isochrones(Z_isocs, m_isocs, Mv_isocs, Te_isocs, nintpl = 10000)
    tf_intpl = time()
    print('Interpolation process took {0} sec.\n'.format(tf_intpl-t0_intpl))

    print('Loading observational data')
    t0_load_data = time()
    Z_obs = [0.02]*26
    Mv = [6, 5.5, 5, 5, 5, 4.5, 4.5, 4.5, 4, 4, 4, 4, 3.5, 3.5, 3.5, 3.5, 2.5, 2.5, 2.5, 2.5, 2.5, 1.5, 1.5, 1.5, 1.5, 1.5]
    Te = [5000,5300, 5450, 5550, 5650, 5600, 5750, 5900, 5000, 5420, 5840, 6250, 4850, 5430, 6010, 6600, 4800, 5435, 6070, 6705, 7350, 4700, 5375, 6050, 6725, 7400]
    obj_id = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    obj_id.reverse()

    Te_err = [100]*26
    Mv_err = [0.1]*26
    Z_err = [0.003]*26
    tf_load_data = time()
    print('Loading data took {0} sec.\n'.format(tf_load_data-t0_load_data))

    print('Calculating G(t)')
    t0_G = time()
    tseq = np.arange(0.1, 13.5, 0.1)
    G = np.zeros((len(Mv), len(tseq)))
    for i in range(len(Mv)):
        print('For star {0}'.format(obj_id[i]))
        G[i] = G_for_loaded_isocs(Te_obs=Te[i], Te_err=Te_err[i],
                                  Mv_obs=Mv[i], Mv_err=Mv_err[i],
                                  Z_obs=Z_obs[i], Z_err=Z_err[i],
                                  Te_isocs=Te_isocs, Mv_isocs=Mv_isocs,
                                  Z_isocs=Z_isocs, m_isocs=m_isocs,
                                  tseq=tseq)
    tf_G = time()
    print('Calculating G(t) for {0} stars took {1} sec.\n'.format(len(Mv), tf_G - t0_G))

    tf_total = time()
    print('Total process took {0} sec.\n'.format(tf_total-t0_total))

    print('Ploting results')
    for i in range(len(Mv)):
        print('For star {0}'.format(obj_id[i]))
        plt.plot(tseq, G[i,:])
        plt.show()

###########
# Testing
#################
#
if 0:
    Te_obs = 4700
    Te_err = 150
    Mv_obs = 1.5
    Mv_err = 0.1
    Z_obs = 0.020
    Z_err = 0.003

    G(Te_obs=Te_obs, Te_err=Te_err, Mv_obs=Mv_obs, Mv_err=Mv_err, Z_obs=Z_obs, Z_err=Z_err,
      path="C:\\Users\\Felipe\\Documents\\Doutorado\\Projeto Idade Isocronal\\Isochrone\\Isochrones_set_parsec1.0_res0.1Ga\\",
      name_base = "parsec1.0_res0.5",
      nintpl = 10000)

if 0:
    t0 = time()
    # Test Stars
    Z_obs = 0.02
    Mv = [6, 5.5, 5, 5, 5, 4.5, 4.5, 4.5, 4, 4, 4, 4, 3.5, 3.5, 3.5, 3.5, 2.5, 2.5, 2.5, 2.5, 2.5, 1.5, 1.5, 1.5, 1.5, 1.5]
    Te = [5000,5300, 5450, 5550, 5650, 5600, 5750, 5900, 5000, 5420, 5840, 6250, 4850, 5430, 6010, 6600, 4800, 5435, 6070, 6705, 7350, 4700, 5375, 6050, 6725, 7400]
    obj_id = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    obj_id.reverse()

    Te_err = 100
    Mv_err = 0.1
    Z_err = 0.003

    N_obj = len(obj_id)
    t_seq = np.arange(0.1, 13.5, 0.1)
    G_all = np.zeros([N_obj, len(t_seq)])

    for i in range(N_obj):
        print 'Obtainning age for the star {0}'.format(obj_id[i])
        print 'Te = {0} | Mv = {1} | Z = {2}'.format(Te[i], Mv[i], Z_obs)

        Te_obs = Te[i]
        Mv_obs = Mv[i]

        G_all[i,] = G(Te_obs=Te_obs, Te_err=Te_err, Mv_obs=Mv_obs, Mv_err=Mv_err, Z_obs=Z_obs, Z_err=Z_err,
                      path="C:\\Users\\Felipe\\Documents\\Doutorado\\Projeto Idade Isocronal\\Isochrone\\Isochrones_set_parsec1.0_res0.1Ga\\",
                      name_base = "parsec1.0_res0.5",
                      nintpl = 500)

    tf = time()
    print 'Total time: {0} sec.'.format(tf-t0)
    np.save(file='Testes_19-02-16.npy', arr=G_all)

