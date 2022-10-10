############################## MAZEPIN v0.2.2 #################################
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
import scipy.integrate   as sci
import os

from   scipy.fftpack     import fft, fftfreq, fftshift
from   mazepin_aux       import *

# I will define a sequence of 10 colors that I can distinguish without many
# problems. If you are not color-blind, you can comment or ignore the following
# lines and go to the definition of functions

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
    
def CerenkovRing(Xmax, ground, theta, UseSlant = True):
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
    
def cos_localtheta(h, theta, RASPASSHeight = 0.):
    ''' Cosine of trajectory's local zenith angle at height h [km]
    
        theta is PrimaryZenAngle in AIRES (deg)
        
        Accepts RASPASS trajectories
    '''
    
    thetarad = theta *np.pi/180.
    RT = 6370.
    
    return np.sqrt(1.-((RT+RASPASSHeight)/(RT+h)*np.sin(thetarad))**2)

def h_IP(RD, RH, theta):
    ''' Returns height of injection point [km] in a RASPASS trajectory
    
        RD: RASPASSDistance [km]
        
        RH: RASSPASSHeight  [km]
        
        theta: PrimaryZenAngle (deg)
    '''
        
    thetarad = theta *np.pi/180.
    RT = 6370.
    
    return np.sqrt((RT+RH)**2+RD**2+2*RD*(RT+RH)*np.cos(thetarad))-RT

def h_RAS(L, RD, RH, theta):
    ''' Returns height in the atmosphere [km] of the point of the shower axis
        that is a distance L [km] away from the IP (converts traversed 
        distance to height)
        
        L: traversed distance from IP [km]
        
        RD: RASPASSDistance [km]
        
        RH: RASPASSHeight [km]
        
        theta: PrimaryZenAngle [deg]
    '''
    RT = 6370.
    
    v = h_IP(RD, RH, theta)
    c = cos_localtheta(v, theta, RASPASSHeight = RH)
    
    return np.sqrt((RT+v)**2+L**2-2*L*(RT+v)*c)-RT

def Xs_to_dist(X, RD, RH, theta, prec = .01):
    ''' Converts slanted depth [g/cm2] to distance covered along shower axis,
        for a RASPASS trajectory. Approximate result (numerical integration)
    
        X: values of slanted depth [km] (numpy.ndarray)
        
        RD: RASPASSDistance [km]
    
        RH: RASPASSHeight [km]
    
        theta: PrimaryZenAngle [deg]
    
        prec: precission for the result
        
        Uses exponential model for atmosphere (for the moment). My plan is to 
        add the Linsley model.
    '''
    
    RT = 6370.
    thetarad = theta *np.pi/180.
    
    def rho(h):
        return 1.23e-3 * np.exp(-h/8.13)
    
    def integrand(h):
        return rho(h)/np.sqrt(1-((RT+RH)/(RT+h)*np.sin(thetarad))**2)
    
    X      = np.sort(X)    # we sort values (just in case)
    X_last = max(X)        # last value, will be the highest
    
    v      = h_IP(RD, RH, theta) # height of injection point
    c      = cos_localtheta(v, theta, RASPASSHeight = RH)
    hmin   = (RT+v)*np.sqrt(1-c*c)-RT # minimum height
    
    hmin   = hmin if hmin >= 0 else 0
    
    x      = 0.
    L      = 0.
    
    dist   = []
    index  = 0
    
    while x < X_last:
        delta_x = sci.quad(integrand, h_RAS(L+prec, RD, RH, theta), \
                           h_RAS(L, RD, RH, theta))[0]*1e5
        # we find matter traversed when we move a distance prec
        x_new   = x + delta_x
        
        if x_new > X[index]: 
            # X[index] is between x and x_new = x+delta_x
            # we use a linear interpolation between the last two values
            
            distance = L + prec/delta_x * (X[index]-x)
            index   += 1
            dist.append(distance)
        
        x  = x_new
        L += prec
        
    return np.array(dist)
    
def table_finder(tab, part, verbose = False):
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
        
        8 : Number and energy of particles at ground vs shower number
        
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
        
    table = tables[part,tab]
    
    if table == 9999:
        raise TypeError('Table is not available in AIRES 19.04.08')
    
    if verbose:
        print('Requested table: '+dict_tab[str(tab)]+', for '+dict_part[str(part)])
        
    return int(table)

def pathfinder(rootdir, tab, part, sim = [''], sep = [''], verbose = False):
    ''' Returns the list of paths inside rootdir corresponding to the requested
        tables (only one type of table per request).
        
        tab (int): type of table. See mazepin_aux.
        
        part (list of int): types of particles. See mazepin_aux
        
        sim (list of str): constraints for file names (all selected files must 
        these strings). Default '', all files with the right extension are kept.
        I am using my rules for naming files (see notes at GitHub repo). 
        
        sep (list of str): Special distinctions. Output will indicate which files
        contain these special strings. Default '', all files with the right 
        extension are kept
        
        Output (with separators sep = [sep1, sep2]):
            
        [[part1, sep1, path1], [part1, sep2, path2], [part2, sep1, path3], ...]
        
        ** sorry for the spaghetti **
    '''
    paths = []
    for subdir, dirs, files in os.walk(rootdir):
        for s in sep:          # loop over distinctions
            for file in files: # loop over files n rootdir
                for p in part: # loop over requested particles
            
                    table_ext = '.t'+str(tables[p,tab]) 
                    # extension of table tab, particle p
                
                    if file.endswith(table_ext) and all([c in file for c in sim]):
                    # if our file has the right extension and all constraints
                        if s == '' or s in file: #we take into account distinctions
                            paths.append([str(p), s, subdir + os.sep + file])
    
    if len(sep) > 1 and len(paths) != len(part) * len(sep):
        raise TypeError('Some tables are not available inside the main directory')

    if verbose:
        print('Requested table: '+dict_tab[str(tab)]+'\n')
        print('--------------------------\n')
        print('Requested particles: \n')
        for p in part:
            print(dict_part[str(p)]+'\n')
        print('--------------------------\n')
        print('Constraints: \n')
        for c in sim:
            print(c+'\n')
        print('--------------------------\n')
        print('Separations: \n')
        for s in sep:
            print(s+'\n')
        print('--------------------------\n')
        for r in paths:
            print(dict_part[r[0]]+', '+r[1]+': \n'+r[2]+'\n')
        
    return paths
            


