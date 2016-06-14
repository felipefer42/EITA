import numpy as np
import os
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

header = ['MODELL', 'MASS', 'AGE', 'LOG_L', 'LOG_TE', 'LOG_R', 'LOG_RAT', 'M_CORE_HE', 'M_CORE_C', 'H_CEN', 'HE_CEN', 
          'C_cen', 'O_cen', 'LX', 'LY', 'LC', 'LNEUTR', 'L_GRAV', 'H_SUP', 'HE_SUP', 'C_SUP', 'N_SUP', 'O_SUP', 
          'PHASE']

path = "/home/felipe/MEGA/Isocpy/Evolutionary_Tracks"

def get_evotrack_filename(Z, Y, M):
    Z_fmt = str(Z)
    Y_fmt = str(Y)
    OUTA = '1.77' if M <= 0.7 else '1.74'
    return "Z{:s}Y{:s}OUTA{:s}_F7_M{:07.3f}.DAT".format(Z_fmt, Y_fmt, OUTA, M)

def get_evotrack_folder(Z, Y):
    Z_fmt = str(Z)
    Y_fmt = str(Y)
    return "Z{:s}Y{:s}".format(Z_fmt,Y_fmt)

# Load test
if 0:
    Z = 0.0001
    Y = round(0.2485+1.78*Z,3)
    M = 1

    evotrack_folder = get_evotrack_folder(Z,Y)
    evotrack_filename = get_evotrack_filename(Z, Y, M)
    load_file = path+'/'+evotrack_folder+'/'+evotrack_filename

    print 'path exists: {}'.format(os.path.isdir(path))
    print 'evolutionary folder exists: {}'.format(os.path.isdir(path+'/'+evotrack_folder))
    print 'file exists: {}'.format(os.path.exists(load_file))
    print load_file

    data = np.loadtxt(load_file, skiprows = 1)
    
    plt.scatter(data[:,4], data[:,5], marker = '+', c = data[:,-1])
    plt.xlim(max(data[:,4]), min(data[:,4]))
    plt.show()

    for i in range(data.shape[1]):
        if i != 2:
            plt.scatter(data[:,2], data[:,i], marker = '+', c = data[:,-1])
            plt.xlabel(header[2])
            plt.ylabel(header[i])
            plt.show()

class Evo_Track(object):
    def __init__(self, Z, M, path = '', load = True):
        self.Z = Z
        self.Mini = M
        self.Y = round(0.2485+1.78*self.Z, 3)

        if load == True:
            evotrack_folder = get_evotrack_folder(self.Z, self.Y)
            evotrack_filename = get_evotrack_filename(self.Z, self.Y, self.Mini)

            load_file = path+'/'+evotrack_folder+'/'+evotrack_filename

            print('loading evolutionary track data from file '+load_file)
            evotrack_data = np.loadtxt(load_file, skiprows = 1, usecols = [0,1,2,3,4,5,23])

            self.model = evotrack_data[:,0]
            self.mass = evotrack_data[:,1]
            self.age = evotrack_data[:,2]
            self.logL = evotrack_data[:,3]
            self.logTe = evotrack_data[:,4]
            self.logR = evotrack_data[:,5]
            self.phase = evotrack_data[:,6]

    def interp(self, age, interp_phases = range(1,16)):
        age = np.array(age)
        # Find first age of each phase
        phase_int = list(self.phase.astype('int'))

        phase_init_age = []
        phase_final_age = []

        for i in interp_phases:
            try:
                # index of the first occurency of i in phase_int
                init_age_index = phase_int.index(i) 
                # index of the last occurency of i in phase_int
                final_age_index = len(phase_int) - phase_int[::-1].index(i) - 1 

            except ValueError:
                phase_init_age.append(False)
                phase_final_age.append(False)
                continue

            phase_init_age.append(self.age[init_age_index])
            if final_age_index == len(self.age)-1:
                phase_final_age.append(self.age[final_age_index])
            else:
                phase_final_age.append(self.age[final_age_index+1])

        # Interpolating for each phase
        mass = np.array([])
        logL = np.array([])
        logTe = np.array([])
        logR = np.array([])
        phase = np.array([])

        for i in range(len(interp_phases)):
            if not phase_init_age[i]: continue

            evotrack_age_filter = (self.age >= phase_init_age[i]) & (self.age <= phase_final_age[i])
            age_filter = (age >= phase_init_age[i]) & (age < phase_final_age[i])
            
            age_step = age[age_filter]

            if len(age_step) == 0: continue

            mass_interp = interp1d(self.age[evotrack_age_filter], self.mass[evotrack_age_filter])
            mass = np.concatenate((mass, mass_interp(age_step)))

            logL_interp = interp1d(self.age[evotrack_age_filter], self.logL[evotrack_age_filter])
            logL = np.concatenate((logL, logL_interp(age_step)))

            logTe_interp = interp1d(self.age[evotrack_age_filter], self.logTe[evotrack_age_filter])
            logTe = np.concatenate((logTe, logTe_interp(age_step)))

            logR_interp = interp1d(self.age[evotrack_age_filter], self.logR[evotrack_age_filter])
            logR = np.concatenate((logR, logR_interp(age_step)))

            phase_interp = interp1d(self.age[evotrack_age_filter], self.phase[evotrack_age_filter])
            phase = np.concatenate((phase, phase_interp(age_step)))

        # Create new interpolated evolutionary track

        evotrack = Evo_Track(Z = self.Z, M = self.Mini, load = False)
        evotrack.mass = mass
        evotrack.age = age
        evotrack.logL = logL
        evotrack.logTe = logTe
        evotrack.logR = logR
        evotrack.phase = phase

        return evotrack

if 1: # Testing age interpolation

    teste = Evo_Track(Z = 0.0001, M = 1, path = '/home/felipe/MEGA/Isocpy/Evolutionary_Tracks')
    ages = np.arange(0.01e9, 5.5e9, 0.001e9)

    teste_interp = teste.interp(age = ages)
    plt.plot(teste.logTe, teste.logL, 'r')
    plt.plot(teste_interp.logTe, teste_interp.logL, 'b')
    plt.xlim(max(teste.logTe), min(teste.logTe))
    plt.show()
