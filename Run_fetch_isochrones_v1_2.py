__author__ = 'Felipe'

# Needs the files fetch_isochrones.py and download_from_url.py in the same directory as this one
# Uses the libraries: os, sys, urllib, mechanize, numpy

# Version 1.2
# Compatibility with fetch_isochrones_v1_2.py

from fetch_isochrones_v1_2 import *
import numpy as np

###################################################################################################################
# Downloading data from http://stev.oapd.inaf.it/cgi-bin/cmd (check there for more information about the options) #
###################################################################################################################

# The output is one file per age, containing the correspondent isochrones for all choosen metallicities

##############
# Parameters #
############################################
# folder where data must be saved
path = "/home/felipe/Documents/Isochrones/Sets/Set4"
# Namebase for the files
filename_base = 'isoc'

t0 = 0.1e9  # Initial age
tf = 13.6e9  # Final age
dt = 0.1e9  # Delta age

z0 = 0.0001  # Initial metallicity
zf = 0.0601  # Final metallicity
dz = 0.001  # Delta metallicity

group_metal = True # If True, different metallicities are grouped in the same file

isoc_kind = 'parsec_CAF09_v1.0'  # Evolutionary Track
photsys_version = 'yang'  # Bolometric correction to 'normal stars'
photsys_file = 'tab_mag_odfnew/tab_mag_ubvrijhk.dat'  # Photometric system file
kind_cspecmag = 'aringer09'  # Bolometric correction to C stars
dust_sourceM = 'nodustM'  # Dust Composition for M stars
dust_sourceC = 'nodustC'  # Dust Composition for C stars
extinction_av = '0.0'  # Interstellar Extinction
imf_file = 'tab_imf/imf_chabrier_lognormal.dat'  # Initial mass function file

############################################

if tf == t0:
    if zf == z0:
        filename_base=filename_base+'_Z{0}_t{1}e9.dat'.format(z0, t0/1e9)
        fetch_isochrones(isoc_kind=isoc_kind,
                         photsys_version=photsys_version,
                         photsys_file=photsys_file,
                         kind_cspecmag=kind_cspecmag,
                         dust_sourceM=dust_sourceM,
                         dust_sourceC=dust_sourceC,
                         extinction_av=extinction_av,
                         imf_file=imf_file,
                         sequence_type=0,
                         isoc_age=t0,
                         isoc_z=z0,
                         path=path,
                         filename=filename_base)
    else:
        if group_metal:
            filename_base=filename_base+'_zi{0}_zf{1}_dz{2}_t{3}.dat'.format(z0,zf,dz,t0/1e9)
            fetch_isochrones(isoc_kind=isoc_kind,
                             photsys_version=photsys_version,
                             photsys_file=photsys_file,
                             kind_cspecmag=kind_cspecmag,
                             dust_sourceM=dust_sourceM,
                             dust_sourceC=dust_sourceC,
                             extinction_av=extinction_av,
                             imf_file=imf_file,
                             sequence_type=2,
                             isoc_age=t0,
                             isoc_z0=z0,
                             isoc_z1=zf,
                             isoc_dz=dz,
                             path=path,
                             filename=filename_base)
        else:
            z = np.arange(z0, zf+dz/10, dz)
            Nfiles = len(z)
            N = 1
            for zi in z:
                print 'Fetching metallicity Z = {0}'.format(zi)
                filename_i = filename_base+'_Z{0}_t{1}e9.dat'.format(zi, t0/1e9)
                fetch_isochrones(isoc_kind=isoc_kind,
                                 photsys_version=photsys_version,
                                 photsys_file=photsys_file,
                                 kind_cspecmag=kind_cspecmag,
                                 dust_sourceM=dust_sourceM,
                                 dust_sourceC=dust_sourceC,
                                 extinction_av=extinction_av,
                                 imf_file=imf_file,
                                 sequence_type=0,
                                 isoc_age=t0,
                                 isoc_z=zi,
                                 path=path,
                                 filename=filename_i)

                print 'Downloaded file {0} of {1}\n'.format(N, Nfiles)
                N += 1
else:
    if group_metal:
        fetch_multiple_isochrones_constant_age(t0=t0,
                                               tf=tf,
                                               dt=dt,
                                               z0=z0,
                                               zf=zf,
                                               dz=dz,
                                               isoc_kind=isoc_kind,
                                               photsys_version=photsys_version,
                                               photsys_file=photsys_file,
                                               kind_cspecmag=kind_cspecmag,
                                               dust_sourceM=dust_sourceM,
                                               dust_sourceC=dust_sourceC,
                                               extinction_av=extinction_av,
                                               imf_file=imf_file,
                                               name_base = filename_base,
                                               path = path)
    else:
        z = np.arange(z0, zf+dz/10, dz)
        t = np.arange(t0, tf+dt/10, dt)
        Nfiles = len(z)*len(t)

        N = 1
        for ti in t:
            print 'Fetching age {0}e9 Ga'.format(ti/1e9)
            for zi in z:
                print 'Fetching metallicity Z = {0}'.format(zi)
                filename_i = filename_base+'_Z{0}_t{1}e9.dat'.format(zi, ti/1e9)
                fetch_isochrones(isoc_kind=isoc_kind,
                                 photsys_version=photsys_version,
                                 photsys_file=photsys_file,
                                 kind_cspecmag=kind_cspecmag,
                                 dust_sourceM=dust_sourceM,
                                 dust_sourceC=dust_sourceC,
                                 extinction_av=extinction_av,
                                 imf_file=imf_file,
                                 sequence_type=0,
                                 isoc_age=ti,
                                 isoc_z=zi,
                                 path=path,
                                 filename=filename_i)

                print 'Downloaded file {0} of {1}\n'.format(N, Nfiles)
                N += 1
