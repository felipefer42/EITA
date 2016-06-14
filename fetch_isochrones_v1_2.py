__author__ = 'Felipe'

from download_from_url import *
from time import sleep
from mechanize import Browser
import warnings
import numpy as np

# View form controls
# br = Browser()
# br.open('http://stev.oapd.inaf.it/cgi-bin/cmd')
# br.select_form(nr = 0)
# for control in br.form.controls: print control

#Version 1.2
#Improved user handling of argument sequence_type
#Included function fetch_multiple_isochrones_constant_metal
#Included function fetch_multiple_isochrones_constant_age

def try_except_wait(function, waiting_time = 60, attempts = 10, msg = 'Attempting to download data', **kwargs):
    from time import sleep
    if attempts:
        try:
            print(msg + ' (Attempts left: ' + str(attempts) + ')')
            function(**kwargs)
        except:
            'Waiting ' + str(waiting_time) + ' seconds for the next attempt.'
            sleep(waiting_time)
            try_except_wait(function, waiting_time, attempts-1, msg, **kwargs)
    else:
        function(kwargs)

def fetch_isochrones(isoc_kind='parsec_CAF09_v1.2S',
                     photsys_version='yang',
                     photsys_file='tab_mag_odfnew/tab_mag_ubvrijhk.dat',
                     kind_cspecmag='aringer09',
                     dust_sourceM='nodustM',
                     dust_sourceC='nodustC',
                     extinction_av='0.0',
                     imf_file='tab_imf/imf_chabrier_lognormal.dat',
                     sequence_type= 'single_isochrone',
                     isoc_age=False,
                     isoc_z =False,
                     isoc_z0=False,
                     isoc_z1=False,
                     isoc_dz=False,
                     isoc_lage0=False,
                     isoc_lage1=False,
                     isoc_dlage=False,
                     path='',
                     filename='Isochrone_teste.dat'):

    #Sequence_type = 'single_isochrone', 'sequence_constant_metallicity', 'sequence_constant_age'
    if sequence_type == 'single_isochrone' or sequence_type == 0: sequence_type = 0
    elif sequence_type == 'constant_metallicity' or sequence_type == 1 : sequence_type = 1
    elif sequence_type == 'constant_age' or sequence_type == 2: sequence_type = 2
    else: raise ValueError("Argument sequence_type must be in ('single_isochrone', 'constant_metallicity', "
                           "'constant_age')")

    warnings.simplefilter('always', UserWarning)

    #Handling bad values given for different sequence types
    if sequence_type == 0:
        if not isoc_age: raise ValueError("For sequence_type == 'single_isochrone', argument isoc_age must be provided")
        if not isoc_z: raise ValueError("For sequence_type == 'single_isochrone', argument isoc_z must be provided")
        if any((isoc_z0, isoc_z1, isoc_dz, isoc_lage0, isoc_lage1, isoc_dlage)):
            warnings.warn("For sequence_type == 'single_isochrone', arguments isoc_z0, isoc_z1, isoc_dz, isoc_lage0, isoc_lage1 and isoc_dlage are not used")

    elif sequence_type == 1:
        if not isoc_z: raise ValueError("For sequence_type == 'constant_metallicity', argument isoc_z must be provided")
        if not isoc_lage0: raise ValueError("For sequence_type == 'constant_metallicity', argument isoc_lage0 must be provided")
        if not isoc_lage1: raise ValueError("For sequence_type == 'constant_metallicity', argument isoc_lage1 must be provided")
        if not isoc_dlage: raise ValueError("For sequence_type == 'constant_metallicity', argument isoc_dlage must be provided")
        if any((isoc_age, isoc_z0, isoc_z1, isoc_dz)):
            warnings.warn("For sequence_type == 'constant_metallicity', arguments isoc_age, isoc_z0, isoc_z1, and isoc_dz are not used")

    elif sequence_type == 2:
        if not isoc_age: raise ValueError("For sequence_type == 'constant_age', argument isoc_age must be provided")
        if not isoc_z0: raise ValueError("For sequence_type == 'constant_age', argument isoc_z0 must be provided")
        if not isoc_z1: raise ValueError("For sequence_type == 'constant_age', argument isoc_z1 must be provided")
        if not isoc_dz: raise ValueError("For sequence_type == 'constant_age', argument isoc_dz must be provided")
        if any((isoc_z, isoc_lage0, isoc_lage1, isoc_dlage)):
            warnings.warn("For sequence_type == 'constant_age', arguments isoc_z, isoc_lage0, isoc_lage1, and isoc_dlage are not used")

    #Error raised when too many isochrones are requested
    if sequence_type == 1:
        N_isoc = len(np.arange(isoc_lage0, isoc_lage1, isoc_dlage))
        if N_isoc > 400:
            raise ValueError("you requested too many isochrones ({0}), maximum allowed is 400.\nTry to increase isoc_dlage or lower the difference between isoc_lage0 and isoc_lage1".fotmat(N_isoc))
    elif sequence_type == 2:
        N_isoc = len(np.arange(isoc_z0, isoc_z1, isoc_dz))
        if N_isoc > 400:
            raise ValueError("you requested too many isochrones ({0}), maximum allowed is 400.\nTry to increase isoc_dz or lower the difference between isoc_z0 and isoc_z1".format(N_isoc))

    #print 'Opening browser'
    br = Browser()
    br.open('http://stev.oapd.inaf.it/cgi-bin/cmd')
    br.select_form(nr = 0)

    #print 'Filling form'
    br.form['isoc_kind'] = [isoc_kind]
    br.form['photsys_version'] = [photsys_version]
    br.form['photsys_file'] = [photsys_file]
    br.form['kind_cspecmag'] = [kind_cspecmag]
    br.form['dust_sourceM'] = [dust_sourceM]
    br.form['dust_sourceC'] = [dust_sourceC]
    br.form['extinction_av'] = (extinction_av)
    br.form['imf_file'] = [imf_file]

    br.find_control("isoc_val").items[sequence_type].selected = True

    if sequence_type == 0:
        br.form['isoc_age'] = str(isoc_age)  # Isochrone age
        br.form['isoc_zeta'] = str(isoc_z)  # Isochrone metallicity

    elif sequence_type == 1:
        br.form['isoc_zeta0'] = str(isoc_z)  # Isochrone metallicity
        br.form['isoc_lage0'] = str(isoc_lage0)  # Isochrone log initial age
        br.form['isoc_lage1'] = str(isoc_lage1) # Isochrone log final age
        br.form['isoc_dlage'] = str(isoc_dlage)  # Isochrone log age step

    elif sequence_type == 2:
        br.form['isoc_age0'] = str(isoc_age)  # Isochrone age
        br.form['isoc_z0'] = str(isoc_z0)  # Isochrone initial metallicity
        br.form['isoc_z1'] = str(isoc_z1)  # Isochrone final metallicity
        br.form['isoc_dz'] = str(isoc_dz)  # Isochrone metallicity step

    #print('Submitting form')
    br.submit()

    #print('Downloading data')
    download_link = list(br.links())[0].absolute_url
    geturl(download_link, path+'/'+filename)
    br.close()
    print('File ' + path+'/'+filename + ' created')

#################################
# Test                          #
#################################
if 0:
    fetch_isochrones(isoc_kind='parsec_CAF09_v1.2S',
                     photsys_version='yang',
                     photsys_file='tab_mag_odfnew/tab_mag_ubvrijhk.dat',
                     kind_cspecmag='aringer09',
                     dust_sourceM='nodustM',
                     dust_sourceC='nodustC',
                     extinction_av='0.0',
                     imf_file='tab_imf/imf_chabrier_lognormal.dat',
                     sequence_type= 'constant_age',
                     isoc_age=4e9,
                     isoc_z =False,
                     isoc_z0=0.0001,
                     isoc_z1=0.03,
                     isoc_dz=0.00001,
                     isoc_lage0=False,
                     isoc_lage1=False,
                     isoc_dlage=False,
                     path="C:\Users\Felipe\Documents\Doutorado\Projeto Idade Isocronal\Isochrone\Testes22-02-16",
                     filename='teste2.dat')

############################################
# fetch_multiple_isochrones_constant_metal #
############################################

def fetch_multiple_isochrones_constant_metal(z0=0.01,
                                             zf=0.03,
                                             dz=0.01,
                                             isoc_lage0=6.6,
                                             isoc_lage1=10.13,
                                             isoc_dlage=0.05,
                                             path='',
                                             name_base='teste',
                                             isoc_kind='parsec_CAF09_v1.2S',
                                             photsys_version='yang',
                                             photsys_file='tab_mag_odfnew/tab_mag_ubvrijhk.dat',
                                             kind_cspecmag='aringer09',
                                             dust_sourceM='nodustM',
                                             dust_sourceC='nodustC',
                                             extinction_av='0.0',
                                             imf_file='tab_imf/imf_chabrier_lognormal.dat'):

    z = np.arange(z0, zf+dz, dz)
    files_N = len(z)

    for i in range(files_N):

        filename_i = name_base+'_lage0{0}_lage1{1}_dlage{2}_Z{3}.dat'.format(isoc_lage0, isoc_lage1, isoc_dlage, z[i])
        print '\nFetching isochrones with metallicity Z = '+str(z[i])+'.'

        fetch_isochrones(isoc_kind=isoc_kind, photsys_version=photsys_version, photsys_file=photsys_file,
                  kind_cspecmag=kind_cspecmag, dust_sourceM=dust_sourceM, dust_sourceC=dust_sourceC,
                  extinction_av=extinction_av, imf_file=imf_file, sequence_type=1, isoc_z=z[i],
                  isoc_lage0=isoc_lage0, isoc_lage1=isoc_lage1, isoc_dlage=isoc_dlage,
                  path=path, filename=filename_i)

        print 'Downloaded file ' + str(i+1) + ' of '+ str(files_N) + '.'

  #################################
  # Test                          #
  #################################
if 0:
    path="C:\Users\Felipe\Documents\Doutorado\Projeto Idade Isocronal\Isochrone\Testes22-02-16"
    fetch_multiple_isochrones_constant_metal(z0=0.01,
                                             zf=0.03,
                                             dz=0.01,
                                             isoc_lage0=6.6,
                                             isoc_lage1=10.13,
                                             isoc_dlage=0.05,
                                             path=path,
                                             name_base='teste',
                                             isoc_kind='parsec_CAF09_v1.2S',
                                             photsys_version='yang',
                                             photsys_file='tab_mag_odfnew/tab_mag_ubvrijhk.dat',
                                             kind_cspecmag='aringer09',
                                             dust_sourceM='nodustM',
                                             dust_sourceC='nodustC',
                                             extinction_av='0.0',
                                             imf_file='tab_imf/imf_chabrier_lognormal.dat')

############################################
# fetch_multiple_isochrones_constant_metal #
############################################

def fetch_multiple_isochrones_constant_age(t0=1e9,
                                           tf=10e9,
                                           dt=1e9,
                                           z0=0.001,
                                           zf=0.03,
                                           dz = 0.01,
                                           isoc_dlage=0.05,
                                           path='',
                                           name_base='teste',
                                           isoc_kind='parsec_CAF09_v1.2S',
                                           photsys_version='yang',
                                           photsys_file='tab_mag_odfnew/tab_mag_ubvrijhk.dat',
                                           kind_cspecmag='aringer09',
                                           dust_sourceM='nodustM',
                                           dust_sourceC='nodustC',
                                           extinction_av='0.0',
                                           imf_file='tab_imf/imf_chabrier_lognormal.dat',
                                           name = 'simple'):

    t = np.arange(t0, tf+dt, dt)
    files_N = len(t)

    for i in range(files_N):

        if name == 'simple':
            filename_i = name_base+'_t{0}.dat'.format(t[i]/1e9)
        else:
            filename_i = name_base+'_z0{0}_zf{1}_dz{2}_t{3}e9.dat'.format(z0, zf, dz, t[i]/1e9)

        print '\nFetching isochrones with age = '+str(t[i]/1e9)+' Ga.'

        fetch_isochrones(isoc_kind=isoc_kind, photsys_version=photsys_version, photsys_file=photsys_file,
                         kind_cspecmag=kind_cspecmag, dust_sourceM=dust_sourceM, dust_sourceC=dust_sourceC,
                         extinction_av=extinction_av, imf_file=imf_file, sequence_type=2, isoc_age=t[i],
                         isoc_z0=z0, isoc_z1=zf, isoc_dz=dz, path=path, filename=filename_i)

        print 'Downloaded file ' + str(i+1) + ' of '+ str(files_N) + '.'

  #################################
  # Test                          #
  #################################
if 0:
    path="C:\Users\Felipe\Documents\Doutorado\Projeto Idade Isocronal\Isochrone\Testes22-02-16"
    fetch_multiple_isochrones_constant_age(t0=1e9,
                                           tf=10e9,
                                           dt=1e9,
                                           z0=0.001,
                                           zf=0.03,
                                           dz = 0.001,
                                           isoc_dlage=0.05,
                                           path=path,
                                           name_base='teste',
                                           isoc_kind='parsec_CAF09_v1.2S',
                                           photsys_version='yang',
                                           photsys_file='tab_mag_odfnew/tab_mag_ubvrijhk.dat',
                                           kind_cspecmag='aringer09',
                                           dust_sourceM='nodustM',
                                           dust_sourceC='nodustC',
                                           extinction_av='0.0',
                                           imf_file='tab_imf/imf_chabrier_lognormal.dat')




def fetch_isochrones_sequence_age_varying_metal(t0 = 0.1e9, ##### Turn obsolete
                                                tf = 15.0e9,
                                                dt = 0.01e9,
                                                z0 = 0.0001,
                                                zf = 0.03,
                                                dz = 0.0001,
                                                isoc_kind='parsec_CAF09_v1.2S',
                                                photsys_version='yang',
                                                photsys_file='tab_mag_odfnew/tab_mag_ubvrijhk.dat',
                                                kind_cspecmag='aringer09',
                                                dust_sourceM='nodustM',
                                                dust_sourceC='nodustC',
                                                extinction_av='0.0',
                                                imf_file='tab_imf/imf_chabrier_lognormal.dat',
                                                filename = [],
                                                path = ''
                                                ):


    number_n = 1
    import time

    # Number of generated files
    files_N = len(np.arange(t0, tf+dt, dt))
    # Dealing with output files names
    if filename == []:
        for i in range(files_N):
            filename.append('isochrone'+str(i))

    if type(filename) == str and files_N > 1:
        filename_temp = []
        for i in range(files_N):
            filename_temp.append(filename+str(i))
        filename = filename_temp

    # Retrieving first file
    print 'Downloading first file to estimate total files size'
    print '>> After this download, your permission will be asked in order to proceed <<'

    print '\nFetching isochrone of '+str(t0/1e9)+' Gy.'

    kwargs = {'isoc_kind': isoc_kind, 'photsys_version': photsys_version, 'photsys_file': photsys_file,
              'kind_cspecmag': kind_cspecmag, 'dust_sourceM': dust_sourceM, 'dust_sourceC': dust_sourceC,
              'extinction_av': extinction_av, 'imf_file': imf_file, 'sequence_type':2, 'isoc_age': str(t0),
              'isoc_z0': str(z0), 'isoc_z1': str(zf), 'isoc_dz': str(dz), 'path':path, 'filename': filename[0]}

    try_except_wait(fetch_isochrones, waiting_time = 60, attempts = 30, msg = 'Attempting to download data',
                    **kwargs)

    print 'Downloaded file number ' + str(number_n) + ' of this section.'
    number_n += 1

    # Checking first file size
    import os
    statinfo = os.stat(path+'\\'+filename[0])
    filesize = statinfo.st_size
    # Estimated total size
    total_size = files_N * filesize

    unit = 'bytes'
    if total_size < 1024:
        pass
    elif (total_size >= 1024 and total_size < 1024**2):
        total_size /= 1024.
        unit = 'Kb'
    elif (total_size >= 1024**2 and total_size < 1024**3):
        total_size /= 1024**2.
        unit = 'Mb'
    else:
        total_size /= 1024**3.
        unit = 'Gb'

    total_size = round(total_size, 2)

    print 'This function will generate approximately ' + str(total_size) + ' ' + unit + ' in output data.'
    proceed = raw_input('Proceed and create the data? (y/n): ')

    if proceed == 'y':
        #print 'Fetching isochrone of '+str(t0/1e9)+' Gy.'
        #print 'File created: ' + path + '//' + filename[0]
        t = np.arange(t0, tf+dt, dt)
        for i in range(1, files_N):
            #if (number_n % 50) == 0 and number_n != 0:
            #    t_wait = 30
            #    time.sleep(t_wait)
            #    print("Waiting " + str(t_wait) + " seconds to next download...")

            print '\nFetching isochrone of '+str(t[i]/1e9)+' Gy.'

            kwargs = {'isoc_kind': isoc_kind, 'photsys_version': photsys_version, 'photsys_file': photsys_file,
                      'kind_cspecmag': kind_cspecmag, 'dust_sourceM': dust_sourceM, 'dust_sourceC': dust_sourceC,
                      'extinction_av': extinction_av, 'imf_file': imf_file, 'sequence_type':2, 'isoc_age': str(t[i]),
                      'isoc_z0': str(z0), 'isoc_z1': str(zf), 'isoc_dz': str(dz), 'path':path, 'filename': filename[i]}

            try_except_wait(fetch_isochrones, waiting_time = 60, attempts = 30, msg = 'Attempting to download data',
                            **kwargs)

            number_n += 1
            print 'Downloaded file number ' + str(number_n) + ' of this section.'

    else:
        print ('Interrupted by the user, no data created')
        os.remove(path+'//'+filename[0])


#Teste
#import numpy as np

#path = 'C:\Users\Felipe\Documents\Doutorado\Isochrone\Isochrones_test1'
#filename = []
#t0 = 6.7e9
#tf = 14.5e9
#dt = 0.1e9
#z0 = 0.0001
#zf = 0.03
#dz = 0.0001

#filename_base = 'test1_t'
#for t in np.arange(t0, tf+dt, dt):
#    filename.append(filename_base + str(round(t/1e9,2)) + 'e9.dat')
#fetch_isochrones_sequence_age_varying_metal(t0 = t0,
#                                            tf = tf,
#                                            dt = dt,
#                                            z0 = z0,
#                                            zf = zf,
#                                            dz = dz,
#                                            isoc_kind='parsec_CAF09_v1.2S',
#                                            photsys_version='yang',
#                                            photsys_file='tab_mag_odfnew/tab_mag_ubvrijhk.dat',
#                                            kind_cspecmag='aringer09',
#                                            dust_sourceM='nodustM',
#                                            dust_sourceC='nodustC',
#                                            extinction_av='0.0',
#                                            imf_file='tab_imf/imf_chabrier_lognormal.dat',
#                                            filename = filename,
#                                            path = path
#                                            )
