############################## MAZEPIN v1.0.0 #################################
''' 
    Welcome to MAZEPIN (Module for an Aires and Zhaires Environment in PythoN)
    
    This module will be developed in order to have an useful library of Python
    functions to work with AIRES and ZHAireS. The current versions (11/2022) 
    of these programs are AIRES 19.04.08 and ZHAireS 1.0.30a (I will also 
    include some utilities to work with RASPASS)
    
    Some of this code will be based on previous programs that I wrote for
    my Master thesis, designed to produce fast graphs and results from
    direct outputs of AIRES and ZHAireS (though it was quite a mess).
    
    The main goal of this module is to have easy-to-use tools in Python to deal
    with the output files, and also to provide the user with a way of producing
    quickly the most relevant plots.
    
    Updates to this module will be uploaded to my PhD_thesis GitHub repository
    (https://github.com/SergioCabana/PhD_thesis) . There you will also find
    a collection of notes with some of the rules that I will be using to name
    tasks in AIRES and ZHAireS
    
    Sergio Cabana Freire. 05/11/2022.
    
    P.S.: Replacing a more reasonable 'Y' with an 'I' in MAZEPIN is 
    absolutely necessary for the joke, though it might go unnoticed if you are 
    not into F1.
'''

# Importing some modules and functions
import numpy             as np
import matplotlib.pyplot as plt
import matplotlib        as mpl
import pandas            as pd
import os

from scipy.fftpack       import fft, fftfreq, fftshift

# I will define a sequence of 10 colors that I can distinguish without many
# problems. If you are not color-blind, you can comment or ignore the following
# lines and go to the definition of functions

colors = ['k', 'royalblue', 'r', 'gold', 'limegreen', 'navy', 'crimson', \
           'turquoise', 'darkorange', 'darkgreen']
    
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=colors)

def gcm2toh(t):
    ''' Converts vertical depth t (g/cm2) to heights in km
        Utilises Linsley parametrization
        h1km, h2km are the inverses of that parametrization
        
        This function was adapted from a C++ implementation by 
        Washington Carvalho
    '''
    
    if t < 0.:
        raise TypeError('Negative atmospheric thickness')
    
    if t < 0.00128293:
        h2km = 1e4*(1.12829e-2-t)
        return h2km
    
    def h1km(t, params):
        a, b, c = params
        return c*np.log(b/(t-a))
    
    params = [ (-1.86556e2, 1.2227e3, 9.9419), (-9.4919e1, 1.1449e3, 8.7815), \
               (6.1289e-1, 1.3056e3, 6.3614),  (0.0, 5.4018e2, 7.7217) ]

    
    if t <= 3.0395:
        return h1km(t, params[3])
    
    elif t <= 271.6991:
        return h1km(t, params[2])
    
    elif t <= 631.1:
        return h1km(t, params[1])
    
    elif t <= 2004.7:
        return h1km(t, params[0])
    
    else:
        raise TypeError('Atmospheric thickness above 2004.647')
        

def htogcm2(h):
    ''' Converts height in km to vertical depth t (g/cm2)
        Utilises Linsley parametrization, t1 y t2
        
        This function was adapted from a C++ implementation by 
        Washington Carvalho
    '''
    
    if h < -5.801:
        raise TypeError('Altitude lower than -5.8 km')
    
    if h > 112.8:
        return 0.
    
    if h >= 100.0:
        t2 = 1.12829e-2-h/1e4
        return t2
    
    def t1(h, params):
        a, b, c = params
        return a+b*np.exp(-h/c)
    
    params = [ (-1.86556e2, 1.2227e3, 9.9419), (-9.4919e1, 1.1449e3, 8.7815), \
               (6.1289e-1, 1.3056e3, 6.3614),  (0.0, 5.4018e2, 7.7217) ]

    
    if h >= 40.0:
        return t1(h, params[3])
    
    elif h >= 10.0:
        return t1(h, params[2])
    
    elif h >= 4.0 :
        return t1(h, params[1])
    
    else:
        return t1(h, params[0])
    
def CherenkovRing(Xmax, ground, theta, UseSlant = True):
    ''' Returns parameters of the Cerenkov elipse at ground level
        Xmax -> Depth of the shower maximum in g/cm2
                UseSlant True  -> Value measured along shower axis
                UseSlant False -> Value measured in the vertical direction
        ground -> Height of ground a.s.l in m
        theta -> Zenith angle in deg
        
        This function was adapted from a C++ implementation by 
        Washington Carvalho
    '''
    R0 = 325.
    
    theta     = theta * np.pi / 180. # conversion to radians
    costet    = np.cos(theta)
    
    Xmax_vert = Xmax 
    
    if UseSlant:
        Xmax_vert = Xmax * costet # vertical depth of shower maximum
    
    alt_max = gcm2toh(Xmax_vert)*1e3
    h_max   = alt_max-ground     # height of shower maximum (vertical)
    d       = h_max/costet       # distance from ground to maximum (along axis)
    
    # refraction index
    Rh = R0 * np.exp( -0.1218 * alt_max / 1e3 )
    nh = 1. + 1e-6 * Rh      # refraction index at the height of the maximum
    
    costetc = 1./nh          # cosine of Cerenkov angle
    
    thetac  = np.arccos(costetc) 
    sintetc = np.sin(thetac)
    
    ########### APPROXIMATE CALCULATION, FIRST OUTPUT #####################
    
    r1 = h_max * ( np.tan(theta+thetac) - np.tan(theta) )
    r2 = h_max * np.tan(thetac) / costet
    
    print('------ INPUTS -------\n')
    if UseSlant:
        print('X_max slanted (g/cm2): %.4f'%Xmax)
    print('X_max vert. (g/cm2): %.4f'%Xmax_vert)
    print('Zenith angle (deg): %.4f'%( theta * 180/np.pi ))
    print('Height of maximum (m): %.4f'%alt_max)
    print('Height of maximum above ground (m): %.4f'%h_max)
    print('Refraction index, nh: %.8f'%nh)
    print('Cerenkov angle (deg): %.4f\n'%( thetac * 180/np.pi ))
    print('--------- APPROXIMATE CALCULATION ----------\n')
    print('r1 (m): %.5f'%r1)
    print('r2 (m): %.5f\n'%r2)
    
    ############ EXACT CALCULATION ######################
    
    cosplus  = np.cos(theta + thetac)
    cosminus = np.cos(theta - thetac)
    a1       = d * sintetc / cosminus
    a2       = d * sintetc / cosplus
    a        = d * np.sin(2.*thetac) * costet / (2.*cosplus*cosminus)
    b        = d * sintetc * costet / np.sqrt(cosplus*cosminus)
    epsilon  = d * sintetc * sintetc * np.sin(theta) / (cosplus*cosminus)
    
    print('------ CALCULO EXACTO ----------\n')
    print('a1 (m): %.4f'%a1)
    print('a2 (m): %.4f'%a2)
    print('a (m): %.4f'%a)
    print('b (m): %.4f'%b)
    print('Eccentricity : %.4f'%epsilon)
    print('D (m): %.4f'%d)
    print('4,5 * r2 (m) = %.4f'%(4.5*r2))
    
def table_finder(tab, part):
    ''' Returns the extension corresponding to each table and particle
    
        tab: type of table (from 0 to 15)
        
        0 : Longitudinal development
        
        1 : Unweighted longitudinal development
        
        2 : Energy longitudinal development
        
        3 : Lateral distribution
        
        4 : Unweighted lateral distribution
        
        5 : Energy distribution at ground
        
        6 : Unweighted energy distribution
        
        7 : Mean arrival time distribution
        
        8 : Number and energy of particles vs shower number
        
        9 : Number of created particles
        
        10: Number of created entries
        
        11: Energy of created particles
        
        12: Longitudinal development of low energy particles
        
        13: Unweighted longitudinal development of low energy particles
        
        14: Energy longitudinal development of low energy particles
        
        15: Longitudinal development of deposited energy
        
        part: type of particle (from 0 to 22)
        
        0 : gammas
        
        1 : electrons
        
        2 : positrons
        
        3 : muons (+)
        
        4 : muons (-)
        
        5 : pions (+)
        
        6 : pions (-)
        
        7 : kaons (+)
        
        8 : kaons (-)
        
        9 : neutrons
        
        10: protons
        
        11: antiprotons
        
        12: nuclei
        
        13: Other charged particles (excluding p and antip)
        
        14: Other neutral particles (excluding neutrons)
        
        15: e+e-
        
        16: mu+mu-
        
        17: pi+pi-
        
        18: K+K-
        
        19: All charged particles
        
        20: All neutral particles
        
        21: All particles
        
        22: All neutrinos
        
        You need the file 'AIRES_table_index.txt'. You can generate it with 
        the script 'Tabe_index_creator.py', uploaded in the GitHub repo
        
        There are some special tables, these are:
        
        0100 : Atmospheric profile
        
        5501 : Xmax and Nmax (charged particles) vs shower number
        
        5511 : First interaction depth and primary energy vs shower number
        
        5513 : Zenith and azimuth vs shower number
    '''
    
    if tab not in range(0,16,1):
        raise TypeError('Please introduce a valid table id between 0 and 15')
    if part not in range(0,23,1):
        raise TypeError('Please introduce a valid particle id between 0 and 22')
        
    table = np.loadtxt('AIRES_table_index.txt')[part,tab]
    
    if table == 9999:
        raise TypeError('Table is not available in AIRES 19.04.08')
        
    return int(table)
    
    
    