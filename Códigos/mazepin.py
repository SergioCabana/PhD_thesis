############################## MAZEPIN v0.9.0 #################################
''' 
    Welcome to MAZEPIN (Module for an Aires and Zhaires Environment in PythoN)
    
    This module will be developed in order to have an useful library of Python
    functions to work with AIRES and ZHAireS. The current versions (11/2022) 
    of these programs are AIRES 19.04.08 and ZHAireS 1.0.30a (I will also 
    include some utilities to work with RASPASS)
    
    Some of this code will be based on previous programs that I wrote for
    my Master thesis, designed to produce fast graphs and results from
    direct outputs of AIRES and ZHAireS (though it was quite a mess).
    
    The main goal of this module is to have tools in Python to deal with the 
    output files, and also to provide the user with a way of producing
    quickly the most relevant plots. Also, I have included a set of functions
    to write input files and preparing simulations to run in remote machines
    with a SGE management system. I would like to include also their adaptation
    to HTCondor.
    
    Updates to this module will be uploaded to my PhD_thesis GitHub repository
    (https://github.com/SergioCabana/PhD_thesis) . There you will also find a
    script with examples on how to use the most important functions.
    
    Disclaimer: This is not a "formal" module, I would rather call it a script
    with a bunch of definitions. To use them fast, the most important thing up to
    now is copy-pasting a huge number of lists with inputs. There are no definitions
    of objects or a well defined structure. Basically, a huge mess that does the work.
    
    Sergio Cabana Freire. 16/11/2022.
    
    P.S.: Replacing a more reasonable 'Y' with an 'I' in MAZEPIN is 
    absolutely necessary for the joke, though it might go unnoticed if you are 
    not into F1.
'''

# Importing some modules and functions
import numpy             as np
import matplotlib.pyplot as plt
import matplotlib        as mpl
import pandas            as pd
# import scipy.integrate   as sci
import os
import paramiko
import time

from   scipy.fftpack     import fft, fftfreq
from   scipy.interpolate import interp1d
from   getpass           import getpass
from   mazepin_aux       import *

# I will define a sequence of 10 colors (mazepin_aux) that I can distinguish 
# without many problems. If you are not color-blind, you can comment or ignore 
# the following line and go to the definition of functions

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=colors)

# RADIUS OF THE EARTH: This is a global definition for the module.
# You can set it to the average value or to the polar value.
# By default, I will use the polar radius

#RT = 6371. # average [km]
RT = 6357. # polar [km]
 
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
    
    return np.sqrt(1.-((RT+RASPASSHeight)/(RT+h)*np.sin(thetarad))**2)

def RAS_Dist(h, theta, RH = 0):
    ''' Returns the appropriate RASPASSDistance for a shower
        with injecton height h [km] and zenith theta [deg]
        In normal cases, RASPASSHeight = 0 km, but you can set it
        to whichever value you need [km]
    '''
    
    c  = np.cos( theta *np.pi/180 )
     
    return -(RT+RH)*c + np.sqrt((RT+RH)**2*c**2 - (RH**2-h**2+2*RT*(RH-h)))

def h_IP(RD, RH, theta):
    ''' Returns height of injection point [km] in a RASPASS trajectory
    
        ==============================================================
        
        RD: RASPASSDistance [km]
        
        RH: RASSPASSHeight  [km]
        
        theta: PrimaryZenAngle (deg)
        
        ==============================================================
    '''
        
    thetarad = theta *np.pi/180.
    
    return np.sqrt((RT+RH)**2+RD**2+2*RD*(RT+RH)*np.cos(thetarad))-RT

def h_RAS(L, RD, RH, theta):
    ''' Returns height in the atmosphere [km] of the point of the shower axis
        that is a distance L [km] away from the IP (converts traversed 
        distance to height)
        
        ====================================================================
        
        L: traversed distance from IP [km]
        
        RD: RASPASSDistance [km]
        
        RH: RASPASSHeight [km]
        
        theta: PrimaryZenAngle [deg]
        
        ====================================================================
    '''
    
    v = h_IP(RD, RH, theta)
    c = cos_localtheta(v, theta, RASPASSHeight = RH)
    
    return np.abs(np.sqrt((RT+v)**2+L**2-2*L*(RT+v)*c)-RT)

def dist_to_Xs(dist, RD, RH, theta, step = .025):
    ''' Converts distance [km] to slanted depth [g/cm2] covered along shower 
        axis, for a RASPASS trajectory. VERY approximate result 
        (rough and easy numerical integration)
    
        ========================================================================
        
        d: values of distance [km] (numpy.ndarray)
        
        RD: RASPASSDistance [km]
    
        RH: RASPASSHeight [km]
    
        theta: PrimaryZenAngle [deg]
        
        step: Step of integration. Default 25 m
        
        ========================================================================
        
    '''    
    
    thetarad = theta *np.pi/180.

    def integrand(h):
        return rho_Lin(h)/np.sqrt(1-((RT+RH)*np.sin(thetarad)/(RT+h))**2)
        # rho_Lin is the Linsley atmosphere. Defined in mazepin_aux
        
    dist    = np.sort(dist)
    
    hmin    = (RT+RH)*np.sin(thetarad) - RT
    hmin    = hmin if hmin > 0 else 0  #minimum height in atmosphere
    
    Xs      = [0.]
    d_start = 0.
    
    for d in dist:
        
        disc = np.arange(d_start, d+step, step) # discretization of interval
        
        x    = Xs[-1]    # last result (we start from here)
        
        for j in range(len(disc)-1):
            h1    = h_RAS(disc[j], RD, RH, theta)
            h2    = h_RAS(disc[j+1], RD, RH, theta)
            
            if hmin < h1 or hmin > h2:
                delta_x = abs(integrand(.5*(h1+h2))*(h2-h1))*1e5
            else:
                delta_x = abs(integrand(.5*(hmin+h1))*(hmin-h1))*1e5 + \
                          abs(integrand(.5*(h2+hmin))*(h2-hmin))*1e5
            
            x += delta_x
        
        Xs.append(x)
        d_start = d
        
    return np.array(Xs[1:])
       

def atmos_size(RH, theta, atmos_height = 100, stop = 'atmos'):
    ''' Returns total distance and traversed matter for a RASPASS trajectory
        starting at the top ot the atmosphere (fixes RASPASSDistance so that 
        happens)
        
        ====================================================================
        
        RH : RASPASSHeight [km]
        
        theta : PrimaryZenAngle [deg]
        
        atmos_height : height of the atmosphere [km] and starting
            point of the trajectory (default 100 km)
            
        stop : point up to which we measure distance and traversed matter
            'atmos' : top of the atmosphere
            'zaxis' : crossing of shower axis with the vertical of observer
            'hmin'  : minimum height reached in the atmosphere
            
        ===================================================================
    '''
    
    RD     = RAS_Dist(atmos_height, theta, RH = RH)  # RASPASSDistance
    
    c      = cos_localtheta(atmos_height, theta, RASPASSHeight = RH)
    Lmin   = (RT+atmos_height) * c    # distance where minimum height is reached
    L_exit = 2*Lmin              # distance where atmos_height is reached again
    
    hmin  = (RT+RH)*np.sin(theta*np.pi/180.)-RT
    
    hmin = hmin if hmin > 0 else 0
    
    if hmin == 0:  # NOT RASPASS !!!!
        Lmin   = RD
        L_exit = RD
        
    if stop == 'atmos':
        return L_exit, dist_to_Xs(np.array([L_exit]), RD, RH, theta)[0]
    elif stop == 'zaxis':
        return RD, dist_to_Xs(np.array([RD]), RD, RH, theta)[0]
    elif stop == 'hmin':
        return Lmin, dist_to_Xs(np.array([Lmin]), RD, RH, theta)[0]
    
    else: 
        raise TypeError(r'Valid stops are: "atmos", "zaxis" and "hmin" ')
        
        
def Xs_to_dist(X, RD, RH, theta, atmos_height = 100):
    ''' Converts slanted depth [g/cm2] to distance covered along shower axis,
        for a RASPASS trajectory. VERY approximate result 
        (rough and easy numerical integration). Inverts previous function.
    
        Interpolation may not be reliable at all, use with caution
        
        ================================================================
        
        X: values of slanted depth [g/cm2] (numpy.ndarray)
        
        RD: RASPASSDistance [km]
    
        RH: RASPASSHeight [km]
    
        theta: PrimaryZenAngle [deg]
        
        =================================================================

    '''
    
    RD_top = RAS_Dist(atmos_height, theta, RH = RH)
    
    offset_dist = RD - RD_top if RD - RD_top > 0 else 0 # distance outside atmosphere

    L_full, X_full = atmos_size(RH, theta)
    
    d = np.linspace(0, L_full+offset_dist, 1000)
    
    Xs = dist_to_Xs(d, RD, RH, theta)
    
    inverse = interp1d(Xs, d)
    
    return np.array([inverse(x) if x < X_full else L_full+offset_dist for x in X])
    
def codes():
    ''' Returns codes (in mazepin) for AIRES tables and particles
    '''
    
    return dict_tab, dict_part

def table_finder(tab, part, verbose = False):
    ''' Returns the extension corresponding to each table and particle
    
        ==============================================================
        
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
            
            16: Energy (per particle) longitudinal development
        
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

        ====================================================================
        
        There are some special tables, these are:
        
            0100 : Atmospheric profile
            
            5501 : Xmax and Nmax (charged particles) vs shower number
            
            5511 : First interaction depth and primary energy vs shower number
            
            5513 : Zenith and azimuth vs shower number
    '''
    
    if tab not in range(0,17,1):
        raise TypeError('Please introduce a valid table id between 0 and 16')
    if part not in range(0,23,1):
        raise TypeError('Please introduce a valid particle id between 0 and 22')
        
    table = [tables[part,0], tables[part, 2]] if tab == 16 else [tables[part, tab]]  # from mazepin_aux
    
    if table == 9999:
        raise TypeError('Table is not available in AIRES 19.04.08')
    
    if verbose:
        print('Requested table: '+dict_tab[str(tab)]+', for '+dict_part[str(part)])
        if tab == 16:
            print('Two tables needed for energy per particle !!')
        # dictionaries at mazepin_aux
        
    return [int(t) for t in table] 

def pathfinder(rootdir, tabs, part, sim = [''], sep = [], out = [], verbose = False):
    ''' Returns the list of paths inside rootdir corresponding to the requested
        tables.
        
        ===================================================================
        
        tabs (list of int): types of tables. See mazepin.table_finder() .
        
        part (list of int): types of particles. See mazepin.table_finder() .
        
        sim (list of str): constraints for file names (all selected files must 
        these strings). Default '', all files with the right extension are kept.
        I am using my rules for naming files (see notes at GitHub repo). 
        
        sep (list of lists of str): Distinctions. Output will indicate which files,
        among the selected with sim, contain these special strings. 
        Default [], all files with the constraint and right extension are kept
        
        out (list of str): Strings that we do not want our files to have
        Might be useful when we have avergaed tables and individual shower tables,
        that will have and extra "s00..." by default
        Default [], no strings are excluded
        
        ===================================================================
        
        Output (with separators sep = [['s1', 's2'], ['s3', 's4']]):
            
        [[tab1, part1, 's1, s2', path1], [tab1, part1, 's3, s4', path2], 
         [tab1, part2, 's1, s2', path3], [tab1, part2, 's3, s4', path4],
         [tab2, part1, 's1, s2', path5], ...]
        
        ** sorry for the spaghetti, I know it is a huge mess :_( **
    '''
    paths    = []
    separate = True if sep != [] else False
    
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:     # loop over files in rootdir
            for tab in tabs:   # loop over requested tables
                for p in part: # loop over requested particles
                
                    if tab != 16:  # NORMAL AIRES TABLES (JUST LOOP OVER FILES AND CHECK EXTENSIONS)
                    
                        table_ext = '.t'+str(tables[p,tab]) 
                        # extension of table tab, particle p (table at mazepin_aux)
                        sim_check = all([c in file for c in sim])
                        out_check = all([c not in file for c in out])
                        
                        if file.endswith(table_ext) and sim_check and out_check:
                        # if our file has the right extension and all constraints
                            if separate:
                                for s in sep: # loop over separations
                                    string = ", ".join([str(item) for item in s])              
                                    if all([item in file for item in s]): #we take into account distinctions
                                        paths.append([str(tab), str(p), string, subdir + os.sep + file])
                                        # save type of table, particle, separation and path
                            else:
                                paths.append([str(tab), str(p), '', subdir + os.sep + file])
                                    
                    elif tab == 16: # energy per particle, we need 2 tables to calculate
                        
                        table_ext1 = '.t'+str(tables[p,0]) 
                        table_ext2 = '.t'+str(tables[p,2])
                        # extension of tables 0 and 2, particle p (table at mazepin_aux)
                        
                        sim_check = all([c in file for c in sim])
                        out_check = all([c not in file for c in out])
                        
                        if file.endswith(table_ext1) and sim_check and out_check:
                        # if our file has the right extension and all constraints
                            if separate:
                                for s in sep: # loop over separations
                                    string = ", ".join([str(item) for item in s])
                                    if all([item in file for item in s]): #we take into account distinctions
                                        filename = file[:-6] 
                                        #name of task, should have both extensions inside rootdir
                                        paths.append([str(tab), str(p), string, \
                                                     [subdir + os.sep + file, subdir + os.sep + filename+table_ext2]])
                                        # save both paths in the same element
                            else:
                                filename = file[:-6] 
                                #name of task, should have both extensions inside rootdir
                                paths.append([str(tab), str(p), '', \
                                             [subdir + os.sep + file, subdir + os.sep + filename+table_ext2]])
                                    
    if len(sep) > 1 and len(paths) < len(part) * len(sep) * len(tabs):
        raise TypeError('Some tables are not available inside the main directory')
    elif len(sep) > 1 and len(paths) > len(part) * len(sep) * len(tabs):
        raise TypeError('Simulations are not constrained enough (add strings to sim or out)')
        
    if separate:
        order = []
        for p in part:
            for s in sep:
                string = ", ".join([str(item) for item in s])
                for elem in paths:
                    if str(p) == elem[1] and string == elem[2]:
                        order.append(elem) 
                        
                        
        paths = order # this is a trick to make the lists in paths ordered as sep is
    
    if verbose:  # all dictionaries in mazepin_aux
        print('Requested tables: \n')
        for t in tabs:
            print(dict_tab[str(t)]+'\n')
        print('--------------------------\n')
        if 16 in tabs:
            print('WARNING: Energy per particle requested. Exporting tables 0 and 2 together\n')
            print('--------------------------\n')
        print('Requested particles: \n')
        for p in part:
            print(dict_part[str(p)]+'\n')
        print('--------------------------\n')
        print('Constraints: \n')
        for c in sim:
            print(c+'\n')
        print('--------------------------\n')
        print('Unwanted strings: \n')
        for s in out:
            print(s+'\n')
        print('--------------------------\n')
        print('Separations: \n')
        for s in sep:
            print(s+'\n')
        print('--------------------------\n')
        for r in paths:
            print(dict_tab[r[0]]+', '+dict_part[r[1]]+', '+r[2]+':\n')
            print(r[3]); print('\n')
        
    return paths


def readfile(filename):
    ''' Reads the AIRES output given. If ground position is specified in 
        the file comments, returns it 
    '''
    file = open(filename, 'r')
    lines = file.readlines()
    file.close()
    
    ground = 'None'
    
    dataset = []
    for line in lines:
        if line[0:5] == '# GRD':
            ground = line.split()
        elif line[0] != '#':
            dataset += [line.split()]
    
    data = np.array(dataset, dtype = float)
    
    ground = float(ground[2]) if ground != 'None' else 'None'
    
    return data[:,1:], ground


def Aires_data(data, error_type = 'sigma', UG = False, slant = False, \
               RASPASS = False, Distance = False, first_int = True, \
               traject = []):
    ''' Returns data to be plotted, and a collection of suitable labels
    
        ===================================================================
        
        data is a single instance of the output of maz.pathfinder, i.e.,
        [tab, part, sep, path]

        error_type (sigma, RMS): which error to return
        
        UG        : Adapts x_axis data for upgoing showers (basically turns it)
        
        slant     : True if we are using slant depth option for export tables
        
        RASPASS   : Adapts x_axis data for RASPASS trajectories
        
        Distance  : Turns x_axis data to distance measured along axis
        
        first_int : Offset in x_axis for different injection heights
        
        traject   : trajectory parameters. [inj_h, theta (0-90)] (not RASPASS), or
                    [RDistance, RHeight, theta] (RASPASS)
        
        ===================================================================
        
        Output:
            
        [xdata, ydata, error, label], x_axis_label, y_axis_label
        
        **WARNING** Table 5001 is not implemented
    '''
    
    tab, part, sep, file = data
    
    EperPart = True if tab == '16' else False # check if we have to find energy per particle
    
    label = dict_part_tex[part] + ', ' + sep
    
    
    ##################### DATA READING ################################
    if not EperPart:
        
        data, grd = readfile(file)
        
        xdata = data[:,0]
        
        ydata = data[:,1]
        
        try:
            if error_type == 'sigma':
                err = data[:,3]
            elif error_type == 'RMS':
                err = data[:,2]
            else:
                raise TypeError('Please introduce a valid error_type ("sigma" or "RMS")')
        
        except IndexError:
            print('Ignore this message if you requested individual shower tables')
            err = np.zeros(np.shape(xdata))
        
    else:
        data1, grd = readfile(file[0]) # longitudinal development
        data2, _  = readfile(file[1])  # energy longitudinal development data
        
        xdata = data1[:,0]; y1 = data1[:, 1]; y2 = data2[:,1]
        ydata = [a/b if b!=0 else 0 for a, b in list(zip(y2, y1))] # divide both datasets
        
        try:
            if error_type == 'sigma':
                err1 = data1[:,3]; err2 = data2[:,3]
                
                rerr1 = [e/v if v !=0 else 0 for e,v in list(zip(err1, y1))]
                rerr2 = [e/v if v !=0 else 0 for e,v in list(zip(err2, y2))]
                
                err = [v * np.sqrt(re1**2+re2**2) for v, re1, re2 in list(zip(ydata, rerr1, rerr2))]
                # error propagation
                
            elif error_type == 'RMS':
                err1 = data1[:,2]; err2 = data2[:,2]
                
                rerr1 = [e/v if v !=0 else 0 for e,v in list(zip(err1, y1))]
                rerr2 = [e/v if v !=0 else 0 for e,v in list(zip(err2, y2))]
                
                err = [v * np.sqrt(re1**2+re2**2) for v, re1, re2 in list(zip(ydata, rerr1, rerr2))]
                # error propagation
            
            else:
                raise TypeError('Please introduce a valid error_type ("sigma" or "RMS")')
        
        except IndexError:
            print('Ignore this message if you requested individual shower tables')
            err = np.zeros(np.shape(xdata))
        
    x_axis_label, y_axis_label = dict_tab_xlab[tab], dict_tab_ylab[tab] #default
    
    ######################### TRAJECTORY PARAMETERS ####################
    traj_input = False
    
    if RASPASS and len(traject) == 3:
        
        RD, RH, theta = traject
        inj_h         = h_IP(RD, RH, theta)
        traj_input    = True
    
    elif not RASPASS and len(traject) == 2:
        
        inj_h, theta = traject
        RD, RH       = RAS_Dist(inj_h, theta, RH = 0), 0
        
        RD_top       = RAS_Dist(100, theta, RH = 0) # RASPASSDistance of an inj ponit at atmosphere limit
        traj_input   = True
        
        
    ######################### X AXIS CONVERSIONS ##########################
    if not RASPASS:  # "normal" upgoing or downgoing showers, RH = 0
        
        if Distance and tab in tabs_x_depth: # we want to convert x_data to distances along axis 
            if not UG:
                x_axis_label = r'Dist. to ground [km]' 
            else:
                x_axis_label = r'Dist. from ground [km]' if not first_int else r'Dist. from first interaction [km]'
            
            if not traj_input:
                raise TypeError('Trajectory parameters are needed')
                
            if not slant:
                raise TypeError('Conversor Xs_to_dist works with X_slanted as input')

            grd_km = Xs_to_dist(np.array([grd]), RD_top, 0, theta) if grd != 'None' else RD_top
            # position of ground from atmosphere limit, along axis
            
            xdata  = Xs_to_dist(xdata, RD_top, 0, theta)
            
            if first_int and UG:
                xdata  =  grd_km - np.array(xdata) - RD
            else:
                xdata = grd_km - np.array(xdata)
                
        elif slant and not Distance and tab in tabs_x_depth: # slant depths in x axis
            
            x_axis_label = r'$X_s$ [$\mathrm{g/cm^2}$]' 
            
            if first_int:
                
                if not traj_input:
                    raise TypeError('Trajectory parameters are needed for first_int')
                    
                start_depth = dist_to_Xs(np.array([RD_top-RD]), RD_top, 0, theta)
                # we start counting at the injection height, upwards or downwards
                xdata = np.array([start_depth-x for x in xdata]) if UG else np.array([x-start_depth for x in xdata]) 
            
            else:
                xdata = np.array([grd-x for x in xdata]) if UG else np.array(xdata)
            

            
    else: # RASPASS SHOWERS
    
        if Distance and tab in tabs_x_depth:
        
            x_axis_label = r'Dist. to Z axis crossing [km]' 
            
            if not traj_input:
                raise TypeError('Trajectory parameters are needed')
                
            if slant:
                raise TypeError('RASPASS tables by default export distance in x axis. Do not use slant')
            
            xdata = np.array([RD-x for x in xdata])
            
        elif slant and not Distance and tab in tabs_x_depth:
            
            x_axis_label = r'$X_s$ [$\mathrm{g/cm^2}$]' 
            

    return [xdata, ydata, err, label], x_axis_label, y_axis_label


def Aires_Plot(input_data, error_type = 'sigma', UG = False, slant = False, \
               RASPASS = False, Distance = False, first_int = True, \
               trajects = [], xlim = [], ylim = [], xscale = 'linear', \
               yscale = 'linear', graph_type = 'step', size = 12, legend = True,\
               title = '', loc_leg = 'best', fmt = '-', marker = 'o', \
               linewidth = 2., alpha = .7, figsize = (7,5)):
    
    
    ''' Plots AIRES outputs, preprocessed using pathfinder and Aires_data
    
        ====================================================================
        
        input_data : output of pathfinder. The selection of graphs is done there
        
        Allows to plot different tables in the same canvas, as long as it is possible
        (i.e., same X_axis coordinate)
        
        error_type, UG, slant, RASPASS, Distance, first_int, trajects
        are Aires_data options
        
        If traject is not needed, place an empty list [] for each plot
        
        xlim, ylim (list): limits for axes
        
        xscale, yscale ("linear", "log"): scales for axes
        
        graph_type: Valid types are "step", "errorbar" or "fill"
        
        size : size of text
        
        legend (bool): sets legend in graph
        
        title (str): title of graph
        
        loc_leg : position of legend
        
        fmt : if graph_style == errorbar, format of plot
        
        marker : if graph_style == errorbar, marker of points
        
        linewidth : for both types of graphs
        
        alpha : opacity of filling color, only for graph_type = "fill"
        
        figsize: typical option of matplotlib.pyplot.figure
            
        =====================================================================
        
        Automatically sets labels and legends in the plot (If legend = True).
        If you want to change labels or remove the legend after calling the4
        function, use returned axes:
            
            axes.legend(['Label for curve 1', 'Label for curve 2', ...])
            axes.get_legend.remove()
            
        For the last operation, you might get AttributeError, just handle the exception:
        
        for a in axes:
            try:
                a.get_legend().remove()
            except AttributeError:
                continue
    '''
    
    tabs = [data_path[0] for data_path in input_data]
    
    x_depth = all([t in tabs_x_depth for t in tabs])
    x_ener  = all([t in tabs_x_ener for t in tabs])
    x_core  = all([t in tabs_x_core for t in tabs]) # lists in mazepin_aux
    
    if not any([x_depth, x_ener, x_core]):
        raise TypeError('All tables must have the same type of x coordinate')
    
    counter    = 0

    xsets      = []
    ysets      = []
    yerrsets   = []
    graph_labs = []
    ylabs      = []
    
    fig = plt.figure(figsize = figsize)

    ax  = fig.add_subplot(111)
    plt.xticks(fontsize = size)
    plt.yticks(fontsize = size)
    
    # scientific notation for offset text

    # form = mpl.ticker.ScalarFormatter(useMathText=True) 
    # form.set_scientific(True) 
    # form.set_powerlimits((-1,1))
    # ax.yaxis.set_major_formatter(form)
    
    for data_path in input_data:
        dataset, x_axis_label, y_axis_label = Aires_data(data_path, error_type,\
                                              UG, slant, RASPASS, Distance, \
                                              first_int, traject = trajects[counter])
        
        xdata, ydata, err, label = dataset
        
        xsets.append(xdata); 
        ysets.append(ydata);
        yerrsets.append(err);
        graph_labs.append(label); 
        ylabs.append(y_axis_label)
    
        counter += 1
        
        
    ylab_diff = list(dict.fromkeys(ylabs)); multi_y = len(ylab_diff) > 1

    y_axis_label = ''
    
    for y_name in ylab_diff:
        y_axis_label = y_axis_label + y_name + ' or '
        
    y_axis_label = y_axis_label[:-4]
    
    counter = 0
    
    for x, y, err, lab in list(zip(xsets, ysets, yerrsets, graph_labs)):
        
        glabel = dict_tab_yleg[tabs[counter]]+', '+lab if multi_y else lab # dict in mazepin_aux
            
        if graph_type == 'step':
            ax.step(x, y, linewidth = linewidth, where = 'mid', label = glabel)

        elif graph_type == 'errorbar':
            ax.errorbar(x, y, err, fmt = fmt, marker=marker, linewidth=linewidth, label = glabel)
        
        elif graph_type == 'fill':
            ax.errorbar(x, y, fmt = fmt, marker = marker, linewidth = linewidth, label = glabel)
            ax.fill_between(x, y-err, y+err, alpha = alpha)

        else:
            raise TypeError('Valid graph styles are "step", "errorbar" or "fill"')
            
        counter += 1
    
    ax.set_xlabel(x_axis_label, size = size)
    ax.set_ylabel(y_axis_label, size = size)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.set_title(title, size = size)
    
    text = ax.yaxis.get_offset_text()
    text.set_size(size)
    
    if not UG and Distance:
        ax.invert_xaxis()
    
    if UG and not Distance and not slant:
        ax.invert_xaxis()

    
    if legend:
        ax.legend(loc = loc_leg, prop = {'size':size})
        
    if len(xlim) == 2:
        ax.set_xlim(xlim[0], xlim[1])
    
    if len(ylim) == 2:
        ax.set_ylim(ylim[0], ylim[1])
        
    plt.show()
    
    return fig, ax
        

def traject_finder(files, RASPASS = False, UG = False):
    ''' Returns trajecory parameters given the name of the file
        **WARNING** It uses my personal conventions for naming files
        
        =================================================================
        
        files: output of pathfinder
        
        ===================================================================
    '''
    trajects = []
    
    for f in files:
        if type(f[3]) == str:
            table_name = f[3].split('\\')[-1] 
        elif type(f[3]) == list:
            table_name = f[3][0].split('\\')[-1] # just if we have asked for table 16
        
        chars = table_name.split('_')
        
        if RASPASS:
            RD = chars[3][:-2]
            RH = chars[4][:-2]
            theta = chars[6][:-3]
            
            trajects.append([float(RD), float(RH), float(theta)])
            
        elif UG:
            inj_h = chars[3][:-2]
            theta = chars[4][:-3]
            
            trajects.append([float(inj_h), 180-float(theta)]) # upgoing axis has theta < 90 in RASPASS
        
        else:
            inj_h = chars[3][:-2]
            theta = chars[4][:-3]
            
            trajects.append([float(inj_h), float(theta)])
    
    return trajects

def ZHAireS_file_mixer(file, rootdir = '', domain = 't'):
    ''' Combines ZHAireS output files in a single one, in the case that
        several runs with different set of antennas were needed
        
        ================================================================
        
        files (list) : List of file paths to datasets that we want to join
        
        rootdir (str) : Main directory, where the new file will be saved. By default takes execution directory
        
        domain (str) : Type of output file ("t" for TimeFresnel or "f" for FreqFresnel)

        ==================================================================
    '''
    
    datasets = []
    n_ant = 0
    
    for i in range(len(file)):
        try:
            data = np.loadtxt(file[i], comments = '#').T

        except ValueError:
            if rootdir == '':
                print('Reading error, try to introduce rootdir and call again')
            else:
                print('Reading error, working on it ...')
                
                name = 'timefresnel-rootnew_'+str(i+1)+'.dat' if domain == 't' else 'freqfresnel-rootnew_'+str(i+1)+'.dat'
                fileout = rootdir + name
                
                f = open(file[i],'r')
                filedata = f.read()
                f.close()
            
                newdata = filedata.replace("-"," -")
            
                f = open(fileout,'w')
                f.write(newdata)
                f.close()
                print('Corrected file: ', fileout)
                
            data = pd.read_csv(fileout, comment = '#', delimiter = '\s+').to_numpy().T
        
        data[1] = data[1]+n_ant
        
        n_ant = np.max(data[1])
        
        datasets.append(data)
        
        
    new_file = 'timefresnel-root_mix.dat' if domain == 't' else 'freqfresnel-root_mix.dat'
    
    output   = rootdir + new_file 

    
    if domain == 'f':
        
        header = ' -----------Freq Fresnel file----------\n'
        header += ' Combined files: \n'
        for arch in file:
            header += ' '+arch+'\n'
        header += ' \n'
        header += ' Columns: \n'
        header += '        1 - Shower Number\n'
        header += '        2 - Antenna Number\n'
        header += '        3 - Antenna X (m)\n'
        header += '        4 - Antenna Y (m)\n'
        header += '        5 - Antenna Z (m)\n'
        header += '        6 - Frequency #\n'
        header += '        7 - Frequency (MHz)\n'
        header += '        8 - |E| (V/M MHz)\n'
        header += '        9 -  Ex  (V/M MHz)\n'
        header += '        10 - Ey  (V/M MHz)\n'
        header += '        11 - Ez  (V/M MHz)\n'
        header += '\n'
        header += ' new shower    1\n'
        header += '\n' 
        
    elif domain == 't':
        
        header = ' -----------Time Fresnel file----------\n'
        header += ' Combined files: \n'
        for arch in file:
            header += ' '+arch+'\n'
        header += ' \n'
        header += 'Columns: \n'
        header += '        1 - Shower Number\n'
        header += '        2 - Antenna Number\n'
        header += '        3 - Antenna X (m)\n'
        header += '        4 - Antenna Y (m)\n'
        header += '        5 - Antenna Z (m)\n'
        header += '        6 - Time (ns)\n'
        header += '        7 - |A| (V/M)\n'
        header += '        8 -  Ax  (V/M)\n'
        header += '        9 -  Ay  (V/M)\n'
        header += '        10 - Az  (V/M)\n'
        header += '        11 -|E| (V/M)\n'
        header += '        12 - Ex  (V/M)\n'
        header += '        13 - Ey  (V/M)\n'
        header += '        14 - Ez  (V/M)\n'
        header += '\n'
        header += ' new shower    1\n'
        header += '\n' 
        
    full_data = np.hstack(tuple(d for d in datasets))
    
    np.savetxt(output, full_data.T, fmt = '%.3e', header = header)
    
    return output

def ZHAireS_setup(file, time_data = True, plot = False, step_antenna = 5):
    ''' Returns:
        
        numpy.ndarray of data
        A list of antenna labels and positions (AIRES coord system)
        A list of indexes, each one indicates where the data of one antenna starts
        A list of frequencies, if it is the case
        
        Plots the position of antennas with their label if Plot = True
        
        ===================================================================
        
        file (str): path to the ZHAireS output file
        
        time_data (bool): True if TimeFresnel output, False if FreqFresnel Output
        
        plot (bool): whether to plot or not the antenna positions
        
        step_antenna (int): Number of antennas to skip when plotting. For example,
        with step_atenna = 5, will plot antennas 5, 10, 15, 20, ...
        
        ===================================================================
        
    '''
    
    data  = np.loadtxt(file, comments = '#').T
    
    antenna_numbers = list(data[1])
    i_ant           = []
    n_ant           = int(max(antenna_numbers))
    
    ant_data        = []
    
    for i in range(n_ant):
        index = antenna_numbers.index(i+1)
        i_ant.append(index)
        ant_data.append([int(antenna_numbers[index]), data[2][index], data[3][index], 100000. - data[4][index]])

    i_ant.append(len(data.T))
    
    if plot:
        
        fig = plt.figure(1)
        ax = fig.add_subplot(111, projection = '3d')
    
        for label, x, y, z in ant_data[::-step_antenna]:
            
            ax.plot(x, y, z, marker = '1', color = 'k', markersize = 12)
            ax.text(x+10, y+10, z+10, str(int(label)))
            ax.set_xlabel('X [m]')
            ax.set_ylabel('Y [m]')
            ax.set_zlabel('Z [m]')

        plt.show()
        
    if not time_data:
        
        freq_list = [data[6][i] for i in range(i_ant[0], i_ant[1])]
        
        return data, ant_data, i_ant, freq_list
    
    else:
        
        return data, ant_data, i_ant

def ZHAireS_data(data, ant_data, i_ant, plot_request, freq = [], time_data = True, freq_request = []):

    ''' Returns dataset and labels suitable for plotting. Only one request per call

        ===========================================================================

        data, ant_data, i_ant, freq: Output of ZHAireS_setup

        plot_request (list): Element of the form [magnitude, variable, [list of requested antenna indices]]

        time_data (bool) : True (False) if TimeFresnel (FreqFresnel) dataset

        freq_request (list): if not time_data, list of frequency compinents to plot if it is the case


    '''

    mag, var, req_antennas = plot_request

    if time_data and [mag, var] not in [[i, j] for i in valid_mag_t for j in valid_vars_t]:
        raise TypeError('Requested plot is not valid')
    if not time_data and [mag, var] not in [[i, j] for i in valid_mag_f for j in valid_vars_f]:
        raise TypeError('Requested plot is not valid')

    if len(req_antennas) == 0:
        raise TypeError('Please request at least one antena for the plot '+mag+' vs '+var)

    var_2d = True if len(var) == 2 else False

    EFluence = mag[:2] == 'EF'

    if EFluence:
        print('Energy Fluence plots are not implemented yet :_)')

    elif time_data and var == 't': # Plot of signal against time
        curves = []
        for a in req_antennas:
                start_index, end_index = i_ant[a-1], i_ant[a]
                x = data[dict_vars_ZHSt[var]][start_index:end_index]
                y = data[dict_vars_ZHSt[mag]][start_index:end_index]
                
                ax, ay, az = [str(c) for c in ant_data[a][1:]]
                curves.append([x, y, '',  'Antenna ('+ax+', '+ay+', '+az+')'])

            
    elif time_data and var == 'f': # FFT of time data

        curves = []
            
        for a in req_antennas:
            x = data[dict_vars_ZHSt['t']][i_ant[a-1]:i_ant[a]]
            y = data[dict_vars_ZHSt[mag]][i_ant[a-1]:i_ant[a]]
            
            N = len(y) # number of signal points
                
            dT = (x[1]-x[0])*1e-9 #time bin in ns
                
            xf = np.abs(fftfreq(N, dT)*1e-6) #frequencies in MHz
                
            yf = np.abs(fft(y))/N # FFT of data

            ax, ay, az = [str(c) for c in ant_data[a][1:]]
            curves.append([xf, yf, '', 'Antenna ('+ax+', '+ay+', '+az+')'])
                

    elif time_data and var not in ['f', 't']: # Plots against coordinates

        if not var_2d: # single coordinate points
            x = []; y = []

            z0   = 100000. if var == 'z' else 0.
            sign = 1. if var == 'z' else -1.
            
            for a in req_antennas:
                start_index, end_index = i_ant[a-1], i_ant[a]
                
                x.append(z0-sign*data[dict_vars_ZHSt[var]][start_index])
                    
                y.append(max(data[dict_vars_ZHSt[mag]][start_index:end_index]))
                    
            curves = [[np.array(x), np.array(y), '', dict_dirs_ZHS[var]]]

        else:
            x = []; y = []; z = []
            varx, vary = var[0], var[1]
            z0   = 100000. if vary == 'z' else 0.
            sign = 1. if var == 'z' else -1.
            
            for a in req_antennas:

                start_index, end_index = i_ant[a-1], i_ant[a]
                
                x.append(data[dict_vars_ZHSt[varx]][start_index])
                    
                y.append(z0-sign*data[dict_vars_ZHSt[vary]][start_index])

                z.append(max(data[dict_vars_ZHSt[mag]][start_index:end_index]))

            curves = [[np.array(x), np.array(y), np.array(z), '']]

    elif not time_data and var == 'f': # FreqFresnel dataset, plot against frequencies

        curves = []
            
        for a in req_antennas:
            start_index, end_index = i_ant[a-1], i_ant[a]
            x = data[dict_vars_ZHSf[var]][start_index:end_index]
            y = data[dict_vars_ZHSf[mag]][start_index:end_index]
            ax, ay, az = [str(c) for c in ant_data[a][1:]]

            curves.append([x, y, '', 'Antenna ('+ax+', '+ay+', '+az+')'])
            

    else: # FreqFresnel dataset and coordinate plot
        
        curves = []

        if len(freq_request) == 0:
            raise TypeError('Please request at least one frequency for the plot '+mag+' vs '+var)

        if not var_2d:
                
            for f in freq_request:
                    
                x = []; y = []
                    
                for a in req_antennas:

                    index = i_ant[a-1] + freq_request.index(float(f))
                        
                    z0   = 100000. if var == 'z' else 0.
                    sign = 1. if var == 'z' else -1.
                    
                    x.append(z0-sign*data[dict_vars_ZHSf[var]][index])
                    
                    y.append(data[dict_vars_ZHSf[mag]][index])
                        
                
                curves.append([np.array(x), np.array(y), '', str(f)+' MHz, '+dict_dirs_ZHS[var]])

        else:
            varx, vary = var[0], var[1]

            for f in freq_request:
                    
                x = []; y = []; z = []
                    
                for a in req_antennas:
                    
                    index = i_ant[a-1] + freq_request.index(float(f))
                        
                    z0   = 100000. if vary == 'z' else 0.
                    sign = 1. if var == 'z' else -1.
                    
                    x.append(data[dict_vars_ZHSf[varx]][index])
                    
                    y.append(z0-sign*data[dict_vars_ZHSf[vary]][index])
    
                    z.append(data[dict_vars_ZHSf[mag]][index])
                        
                    
                curves.append([x, y, z, str(f)+' MHz'])

    return [mag, var, curves]


def ZHAireS_plot(data, ant_data, i_ant, plots, freqs = [], freq_request = [], \
                 mix_plot = True, time_data = True, lims_cmap = [], fmt = '-',\
                 marker = 'o', linewidth = 2., size = 12, legend = True, \
                 loc_leg = 'best', figsize = (7, 5)):
    
    ''' Represents requested plots
    
        =============================================================================================
        
        data, ant_data, i_ant, freqs: output of ZHAireS_setup
        
        plots (list) : list of elements like [mag, var, [list of antenna indices]]
        
            Valid magnitudes are: 'A', 'Ax', 'Ay', 'Az', 'E', 'Ex', 'Ey', 'Ez', 'EF', 'EFx', 'EFy', 'EFz'
            Valid variables are : 't', 'f', 'x', 'y', 'z', 'xy', 'xz', 'yz'
        
        mix_plot (bool) : If False, each requested plot goes into a different canvas
        
        time_data (bool) : True (False) if we are using TimeFresnel (FreqFresnel) outputs
        
        freq_request (list) : If time_data = False, list of frequency components to plot
        
        lims_cmap (list) : lower and upper limit for cmap (only for plots against 2 coordinates)
        
        fmt, marker, linewidth : usual specifications for plotting
        
        size : fontsize
        
        legend (bool) : Whether to set or not a legend
        
        loc_leg: placing of legend, typical matplotlib options
        
        figsize (tuple) : size of canvas 
        ============================================================================================
        
        Automatically sets labels and legends in the plot (If legend = True).
        If you want to change labels or remove the legend after calling the4
        function, use returned axes:
            
            axes.legend(['Label for curve 1', 'Label for curve 2', ...])
            axes.get_legend.remove()
            
        For the last operation, you might get AttributeError, just handle the exception:
        
        for a in axes:
            try:
                a.get_legend().remove()
            except AttributeError:
                continue
            
        Energy fluence plots are not implemented.
        I would like to include also plots of the Fourier transformed time signal against
        coordinates, for a continuous range of frequency
        
    '''

    n_plots = len(plots)
    
    axes = []
    
    dict_vars = dict_labs_ZHSt_tex if time_data else dict_labs_ZHSf_tex
    dict_legs = dict_labs_ZHSt_tex_2 if time_data else dict_labs_ZHSf_tex_2 
    
    def launch_plot(fig, ax, mag, multi_y, curve, var_2d = False, 
                    lims_cmap = [], fmt = '-', marker = 'o', linewidth = 2.):
        
        ''' Assistant plotting function. Given the figure, axes, magnitude and curve
            will plot with the specified properties.
            Multi_y is a bool to state whether we have more than one variable in the 
            y axis
        '''
        
        x, y, z, label = curve
        
        if var_2d:
            
            df = pd.DataFrame.from_dict(np.array([x, y, z]).T)
            
            df.columns = ['X', 'Y', 'Z']
            
            pivotted = df.pivot('Y', 'X', 'Z').to_numpy()
            
            ny, nx = pivotted.shape
            
            x_min, x_max, y_min, y_max =  min(x), max(x), min(y), max(y)
            
            deltax, deltay = (x_max-x_min)/(nx), (y_max-y_min)/(ny)
            
            extent = [x_min-deltax/2, x_max+deltax/2, y_min-deltay/2, y_max+deltay/2]
            
            if len(lims_cmap) > 0:
                sc = ax.imshow(pivotted , origin = 'lower', cmap = 'gnuplot', 
                               aspect = 'auto', extent = extent, vmin = lims_cmap[0], vmax = lims_cmap[1])
            else:
                sc = ax.imshow(pivotted , origin = 'lower', cmap = 'gnuplot', 
                               aspect = 'auto', extent = extent)
                
            cbar = fig.colorbar(sc)
            cbar.set_label(dict_vars[mag], size = size, rotation = 270, labelpad = 30)
            cbar.ax.tick_params(labelsize=size)
            
            ax.xaxis.set_tick_params(labelsize=size)
            ax.yaxis.set_tick_params(labelsize=size)
            
            if label != '':
                ax.set_title(label, size = size)
        
        elif label != '':
            label = dict_legs[mag] + ', ' + label if multi_y else label
            ax.errorbar(x, y, label = label, fmt = fmt, marker = marker, linewidth = linewidth)
            ax.xaxis.set_tick_params(labelsize=size)
            ax.yaxis.set_tick_params(labelsize=size)
            
        elif multi_y:
            label = dict_vars[mag]
            ax.errorbar(x, y, label = label, fmt = fmt, marker = marker, linewidth = linewidth)
            ax.xaxis.set_tick_params(labelsize=size)
            ax.yaxis.set_tick_params(labelsize=size)
            
        else:
            ax.errorbar(x, y, fmt = fmt, marker = marker, linewidth = linewidth)
            ax.xaxis.set_tick_params(labelsize=size)
            ax.yaxis.set_tick_params(labelsize=size)
            
            
    if mix_plot:
        variables = [ p[1] for p in plots ]
        index_t   = [ i for i in range(n_plots) if plots[i][1] == 't']
        index_f   = [ i for i in range(n_plots) if plots[i][1] == 'f']
        index_c   = [ i for i in range(n_plots) if plots[i][1] in ['x', 'y', 'z']]
        index_2   = [ i for i in range(n_plots) if i not in index_t+index_f+index_c]
        
        if len(index_t) > 0:
            figt = plt.figure(figsize = figsize)
            axt  = figt.add_subplot(111)
            axt.set_xlabel(r'$t$ [$\mathrm{ns}$]', size = size)
            
            variables = [dict_vars[m[0]] for m in [plots[i] for i in index_t]]
            ylab_diff = list(dict.fromkeys(variables)); multi_yt = len(ylab_diff) > 1

            y_axis_label = ''
            for y_name in ylab_diff:
                y_axis_label = y_axis_label + y_name + ' or '
                
            y_axis_label = y_axis_label[:-4]
            
            axt.set_ylabel(y_axis_label, size = size)
            
            axes.append(axt)
            
        if len(index_f) > 0:
            figf = plt.figure(figsize = figsize)
            axf  = figf.add_subplot(111)
            axf.set_xlabel(r'Frequency [$\mathrm{MHz}$]', size = size)
            
            variables = [dict_vars[m[0]] for m in [plots[i] for i in index_f]]
            ylab_diff = list(dict.fromkeys(variables)); multi_yf = len(ylab_diff) > 1

            y_axis_label = ''
            for y_name in ylab_diff:
                y_axis_label = y_axis_label + y_name + ' or '
                
            y_axis_label = y_axis_label[:-4]
            
            axf.set_ylabel(y_axis_label, size = size)
            
            axes.append(axf)
            
        if len(index_c) > 0:
            figc = plt.figure(figsize = figsize)
            axc  = figc.add_subplot(111)
            axc.set_xlabel(r'Distance to core [$\mathrm{m}$]', size = size)
            
            variables = [dict_vars[m[0]] for m in [plots[i] for i in index_c]]
            ylab_diff = list(dict.fromkeys(variables)); multi_yc = len(ylab_diff) > 1

            y_axis_label = ''
            for y_name in ylab_diff:
                y_axis_label = y_axis_label + y_name + ' or '
                
            y_axis_label = y_axis_label[:-4]
            
            axc.set_ylabel(y_axis_label, size = size)
            
            axes.append(axc)
        
        for i in index_t:
            mag, var, curves = ZHAireS_data(data, ant_data, i_ant, plot_request = plots[i])
            
            for c in curves:
                launch_plot(figt, axt, mag, multi_yt, c, fmt = fmt, 
                            marker = marker, linewidth = linewidth)
            
            if legend:
                axt.legend(loc = loc_leg, prop={'size':size})
                
        for i in index_f:
            if time_data:
                mag, var, curves = ZHAireS_data(data, ant_data, i_ant, plot_request = plots[i])
            else:
                mag, var, curves = ZHAireS_data(data, ant_data, i_ant, plot_request = plots[i], freq = freqs,
                                               freq_request=freq_request, time_data = False)
            for c in curves:
                launch_plot(figf, axf, mag, multi_yf, c, fmt = fmt, marker = marker, linewidth = linewidth)
                
            if legend:
                axf.legend(loc = loc_leg, prop={'size':size})
                
        for i in index_c:
            if time_data:
                mag, var, curves = ZHAireS_data(data, ant_data, i_ant, plot_request = plots[i])
            else:
                mag, var, curves = ZHAireS_data(data, ant_data, i_ant, plot_request = plots[i], freq = freqs,
                                               freq_request=freq_request, time_data = False)
            for c in curves:
                launch_plot(figc, axc, mag, multi_yc, c, fmt = fmt, marker = marker, linewidth = linewidth)
                
            if legend:
                axc.legend(loc = loc_leg, prop={'size':size})
                
        for i in index_2:
            
            if time_data:
                mag, var, curves = ZHAireS_data(data, ant_data, i_ant, plot_request = plots[i])
            else:
                mag, var, curves = ZHAireS_data(data, ant_data, i_ant, plot_request = plots[i], freq = freqs,
                                               freq_request=freq_request, time_data = False)
                
            for c in curves:
                fig = plt.figure(figsize = figsize)
                ax  = fig.add_subplot(111)
                ax.set_xlabel(dict_vars[var[0]], size = size)
                ax.set_ylabel(dict_vars[var[1]], size = size)
                axes.append(ax)
                launch_plot(fig, ax, mag, False, c, var_2d = True, lims_cmap = lims_cmap)
            
    
    else:
        for plot in plots:
            mag, var, curves = ZHAireS_data(data, ant_data, i_ant, plot_request = plot, freq = freqs,
                                            time_data = time_data, freq_request = freq_request)
            var_2d = len(var) > 1
            
            if not var_2d:
                fig = plt.figure(figsize = figsize)
                ax = fig.add_subplot(111)
                ax.xaxis.set_tick_params(labelsize=size)
                ax.yaxis.set_tick_params(labelsize=size)
                
                axes.append(ax)
                
                for c in curves:
                    launch_plot(fig, ax, mag, False, c, var_2d, \
                                lims_cmap = lims_cmap, fmt = fmt, marker = marker, linewidth = linewidth)
                        
                ax.set_xlabel(dict_vars[var], size = size)
                ax.set_ylabel(dict_vars[mag], size = size)
                        
            else:
                for c in curves:
                    fig = plt.figure(figsize = figsize)
                    ax = fig.add_subplot(111)
                    ax.xaxis.set_tick_params(labelsize=size)
                    ax.yaxis.set_tick_params(labelsize=size)
                    ax.set_xlabel(dict_vars[var[0]], size = size)
                    ax.set_ylabel(dict_vars[var[1]], size = size)
                    axes.append(ax)
                    launch_plot(fig, ax, mag, False, c, var_2d, \
                                lims_cmap = lims_cmap, fmt = fmt, marker = marker, linewidth = linewidth)
                    
                

            if legend:
                ax.legend(loc = loc_leg, prop={'size':size})
                
    return axes
            

def input_file(task_name, basic, trajectory,  sim_control, RASPASS = False, upgoing = False, \
               ZHAireS = False, ZHAireS_control = [], exports = [], extra = [], set_date = False,\
               save_path = ''):
    
    ''' Creates an AIRES/ZHAireS input file
    
        ==============================================================================
        
        task_name (str): name of task (without .inp extension)
        
        basic (list of str): list of key simulation parameters
        
            [primary, energy (in eV), GeomagneticField, PropagatePrimary, Site, Ground]
            
        trajectory (list of str): depend on primary
            
            [inj_h (in km), theta (in deg), phi (in deg)] for normal AIRES
            
            [ALTITUDE (in km), theta (in deg), phi (in deg)] for upgoing
            
            [RASPASSDistance (in m), RASPASSHeight (in m), RASPASSTimeShift (in ns), theta (in deg), phi (in deg)] for RASPASS
            
        sim_control (list of str): parameters to control the simulation (default options in mazepin_aux):
            
        [TotalShowers, RunsPerProcess, ShowerPerRun, RandomSeed, ObsLevels, ThinningEnergy, ThinnningWFactor, \
         SaveInFile's, electron cuts, gamma cuts, PerShowerData, ExportPerShower]
         
        RASPASS, upgoing (bool): Control of special primaries
        
        ZHAireS: Calculation of radio emission on or off
        
        ZHAireS_control (list of str): ZHAireS instructions:
            [FresnelTime (On/Off), FresnelFreq (On/Off), list_of_antenna_coord (list of [x y z] in m), list_of_freq (list of numbers) in MHz, direct_instructs (list of str)]
            
        exports (list): tables to export. Format is [tab, part, options (str)] with mazepin codes (maz.codes())
        
        set_date (bool): If True, no Date will be specified. This function uses 2010 7 11 as default date,
            but if set_date = False then no directive will be included (unless given in extras)
        save_path: directory to save input file
        
        ============================================================================
    '''
    if len(task_name) > 64:
        raise TypeError('Task Name will get truncated (max. 64 chars)')
        
    if len(basic) != 6:
        raise TypeError('Wrong length of "basic"')

    if RASPASS and len(trajectory) != 5:
        raise TypeError('Wrong length of "trajectory"')
        
    if not RASPASS and len(trajectory) != 3:
        raise TypeError('Wrong length of "trajectory"')
        
    path = os.path.join(save_path, task_name+'.inp')
    
    with open(path, 'w') as f:
        f.write('########## CREATED WITH mazepin.input_file ########### \n')
        f.write('TaskName '+task_name+' \n')
        
        if RASPASS:
            f.write('########### RASPASS PRIMARY ############## \n')
            f.write('AddSpecialParticle RASPASS'+basic[0]+ ' ./RASPASSprimary '+basic[0]+' \n')
            
            basic[0] = 'RASPASS'+basic[0]
            
        elif upgoing:
            f.write('########### UPGOING PRIMARY ############## \n')
            f.write('AddSpecialParticle U'+basic[0]+ ' ./uprimary '+basic[0]+' \n')
            
            basic[0] = 'U'+basic[0]
            
        f.write('########### BASIC: PRIMARY TYPE, ENERGY AND TRAJECTORY ############## \n')
        for i in range(len(basic)):
            if basic[i] != '':    
                f.write(dict_basic_IDL[str(i)] + ' '+ basic[i]+basic_units[str(i)]+' \n')
        
        if not set_date:
            f.write('Date 2010 7 11 \n')
        
        #trajectory
        
        if RASPASS:
            RD, RH, RTs, theta, phi = trajectory
            f.write('SetGlobal RASPASSDistance '+RD+' \n')
            f.write('SetGlobal RASPASSHeight '+RH+' \n')
            f.write('SetGlobal RASPASSTimeShift '+RTs+' \n')
            f.write('PrimaryZenAngle '+theta+' deg \n')
            f.write('PrimaryAzimAngle '+phi+' deg \n')
        
        elif upgoing:
            ALT, theta, phi = trajectory
            f.write('SetGlobal ALTITUDE '+ALT+' \n')
            f.write('PrimaryZenAngle '+theta+' deg \n')
            f.write('PrimaryAzimAngle '+phi+' deg \n')
        
        else:
            inj_h, theta, phi = trajectory
            f.write('InjectionHeight '+inj_h+' km \n')
            f.write('PrimaryZenAngle '+theta+' deg \n')
            f.write('PrimaryAzimAngle '+phi+' deg \n')
            
        f.write('########### SIMULATION CONTROL ############## \n')
        for i in range(len(sim_control)):
            if sim_control[i] != '':    
                f.write(dict_control_IDL[str(i)] + ' '+ sim_control[i]+' \n')
        
        f.write('########### EXTRA INSTRUCTIONS ############## \n')
        for e in extra: 
            f.write(e+' \n')
            
        if ZHAireS:
            
            Time, Freq, Antenas, Freq_list, Direct = ZHAireS_control
            
            f.write('########### ZHAireS CALCULATIONS ############## \n')
            f.write('ZHAireS On \n')
            f.write('FresnelTime '+Time+' \n')
            f.write('FresnelFreq '+Freq+' \n')
            
            i = 0
            for x, y, z in Antenas:
                f.write('AddAntenna '+str(i)+' '+str(x)+ ' ' +str(y) + ' '+str(z)+' \n')
                i+=1
            for fr in Freq_list:
                f.write('AddFrequency '+str(fr)+' \n')
                i+=1
            
            for instruction in Direct:
                f.write(instruction + ' \n')
         
        f.write('########### EXPORTING DATA ############## \n')
        for tab, part, option in exports:

            table = table_finder(tab, part)
            
            if table == 9999:
                raise TypeError('Requested table does not exist')
            
            opt = ' Opt '+option if option != '' else ''
            f.write('ExportTables '+str(table[0])+opt+' \n')
            
        f.write('End')
        f.close()
        
        return path, task_name+'.inp'


def shell_script_SGE(job_name, input_file_dir, input_files, local_dir = '', \
                 main = '/home2/', user = 'sergio.cabana/', \
                 exe = ['aires/bin/ZHAireSRASPASS', 'SpecialPrimaries/RASPASSprimary'],\
                 program = 'ZHAireSRASPASS', autorename = True):
    
    ''' Creates a shell script (.sh extension) that can be submitted to SGE queue 
        systems (at least inside IGFAE nodes, with scratch). 
        It is adapted for the cluster I use but easy to modify
        
        =====================================================================
        
        job_name: job ID for queue system
        
        local_dir: directory to save the created .sh file
        
        input_file_dir: directory of input files (inside main/user/). Ends with /
        
        input_files : list of names of input files
        
        main : default root directory
        
        user : directory of current user
        
        exe  : paths (inside main) to executable binaries to import and use

        program : name of binary executable to run the input files
        
        autorename: bool to rename exported tables with opt a. Use carefully,
            only works for two different export options, namely default and Opt a
            (or Opt p and Opt ap for RASPASS)
        
        =======================================================================
    '''

    path = os.path.join(local_dir, 'script_shell_'+job_name+'.sh')
    
    header = (        
'''############ Created with mazepin.shell_script() ################
#!/bin/bash
#$ -N '''+job_name+'''
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q auger.q
#if [ `hostname -s` -eq `nodo036.inv.usc.es` ]; then
#   sleep 1h
#fi
'''
    )
    
    scratch_dir = '/scratch/' + user + 'Aires_tmp/' + job_name + '/'

    with open(path, 'w') as f:
        
        f.write(header)
        f.write('[[ -d ' + scratch_dir + ' ]] || mkdir -p ' + scratch_dir + ' \n')
        f.write('cd ' + scratch_dir + ' \n')
        
        f.write('############ EXECUTABLES ########### \n')
        for e in exe:
            
            f.write('cp ' + main + user + e + ' ' + scratch_dir + ' \n')
        
        
        f.write('############ INPUT FILES ########### \n')
        for inp in input_files:
            
            f.write('cp ' + main + user + input_file_dir + inp + ' ' + scratch_dir + ' \n')
            
        
        f.write('############ EXECUTION ########### \n')
        for inp in input_files:
            
            f.write('./' + program + ' < ' + inp + ' > ' + inp[:-4] + '.out' + ' \n')
            
        if autorename:
            f.write('############ RENAMING EXTENSIONS ########### \n')
            f.write(rename1 + ' \n')
            f.write(rename2 + ' \n') #loops defined at mazepin_aux
            
        f.write('############ OUTPUT SAVING AND EXIT ########### \n')

        f.write('cp *.* ' + main + user + input_file_dir + ' \n')
        f.write('rm -rf ' + scratch_dir + ' \n')
                    
        f.close()
        
        return path, 'script_shell_'+job_name+'.sh'
    

def setup_simulation_SGE(task_names, basics, trajects, sim_controls, exports, extras, jobID, \
              RASPASS = False, upgoing = False, ZHAireS = False, \
              ZHAireS_control = [], remote_main = '/home2/', user = 'sergio.cabana/', remote_dir = '', \
              exe = ['aires/bin/ZHAireSRASPASS', 'SpecialPrimaries/RASPASSprimary', 'SpecialPrimaries/uprimary'], \
              program = 'ZHAireSRASPASS', autorename = True, local_savepath = '', \
              eco = False, server = 'mastercr1.igfae.usc.es', node = 'nodo014', \
              username = 'sergio.cabana'):
    
    ''' Creates all the files (.inp, .sh) needed to run simulations inside a SGE system,
        and uploads them to your remote machine
        
        ============================================================================
        
        task_names (list): names of tasks (max 64 char), one per each input file required
        
        basics (list): set of basic parameters (arg of mazepin.input_file()). one per each input file
        
        trajects (list): set of trajectory params (arg of mazepin.input_file()). one per each input file
        
        sim_control (list): set of parameters to control simulation (arg of mazepin.input_file()). Default options at mazepin_aux. one per each input file
        
        exports (list): lists of [tab, part, option] to export in each input (arg of mazepin.input_file()). one per each input file.
        
        extras (list): lists of extra instructions to include in each input file (arg of mazepin.input_file()). one per each input file.
        
        jobID (str): Identifier of remote job (qsub)
        
        RASPASS, upgoing (bool): Special particle options
        
        ZHAireS (bool): input files include radio emission calculations
        
        ZHAireS_control: (list): set of parameters to run ZHAireS. (arg of mazepin.input_file()).
        
        remote_main (str, ends with /): root directory in remote machine
        
        user (str, edns with /): remote machine user. All directories in remote machine should be like /remote_dir/user/something
        
        remote_dir (str, ends with /): directory inside remote_dir/user/ to upload files and store outputs
            
        exe (list of str): list of paths inside remote_dir/user/ where required executables are
        
        program: name of executable that runs the simulations
        
        autorename: bool to include loops in shell script to rename tables with Opt a

        local_savepath (str): directory to save created files
        
        eco (bool): True to send only one shell script to remote (otherwise, one per input file to parallelize)
        
        server: host remote machine
        
        node : where to ssh to from remote server
        
        username: name of remote user
        
        ===========================================================================================
        
        Returns the list of paths to uploaded shell scripts inside the remote machine 
        
    '''
    lengths = [len(task_names), len(basics), len(trajects), len(sim_controls), len(exports), len(extras)]
    
    if not all([l == lengths[0] for l in lengths]):
        raise TypeError('Mismatch of lengths of task_names, basics, sim_controls, exports, extras')
        
    if ZHAireS and len(ZHAireS_control) != len(task_names):
        raise TypeError('Length of ZHAireS_control does not match number of input files')
        
    
    # We create all the input files
    simulations = list(zip(task_names, basics, trajects, sim_controls, exports, extras))
    counter     = 0
    
    local_paths = []
    input_files = []
    
    for taskname, basic, traj, sim_control, exp, ext in simulations:
        
        ZHScon = ZHAireS_control[counter] if ZHAireS else []
        
        path, name = input_file(task_name = taskname, basic = basic, trajectory = traj, \
                                sim_control = sim_control, RASPASS = RASPASS, upgoing = upgoing, \
                                ZHAireS = ZHAireS, ZHAireS_control = ZHScon, exports = exp, \
                                extra = ext, save_path = local_savepath)
            
        local_paths.append(path)
        input_files.append(name)
        
        counter += 1
        
    # we create all the required shell_scripts
    
    local_paths_shell = []
    shell_scripts     = []
    
    if eco:
        path, name = shell_script_SGE(job_name = jobID, input_file_dir = remote_dir, \
                                  input_files = input_files, local_dir = local_savepath, \
                                  main = remote_main, user = user, exe = exe, program = program,\
                                  autorename = autorename)
        
        local_paths_shell.append(path)
        shell_scripts.append(name)
            
    else:
        for i in range(len(input_files)):
            path, name = shell_script_SGE(job_name = jobID+str(i), input_file_dir = remote_dir, \
                                      input_files = [input_files[i]], local_dir = local_savepath, \
                                      main = remote_main, user = user, exe = exe, program = program, \
                                      autorename = autorename)
            
            local_paths_shell.append(path)
            shell_scripts.append(name) 
            
    # we connect to our remote host
    
    print('PASSWORD INPUT: Be careful, it might be visible in some terminals')
    print('If you did not get a warning before this, it should work properly')
    
    password = getpass('Please introduce your remote machine password: ')
    
    host = paramiko.SSHClient()
    host.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    host.connect(server, username=username, password=password)

    # we open a communication channel
    hosttransport = host.get_transport()
    dest_addr = (node, 22)
    local_addr = (server, 22)
    hostchannel = hosttransport.open_channel("direct-tcpip", dest_addr, local_addr)

    # connection with node 
    nodehost = paramiko.SSHClient()
    nodehost.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    nodehost.connect(node, username=username, password=password, sock=hostchannel)

    # open sftp channel to upload files
    
    ftp_client  = host.open_sftp()
    
    remote_path = remote_main + user + remote_dir
    
    input_remote_paths = []
    
    shell_remote_paths = []
    
    for inppath, inpname in list(zip(local_paths, input_files)):
        
        ftp_client.put(localpath = inppath, remotepath = remote_path + inpname, confirm = False)
        
        input_remote_paths.append(remote_path + inpname)
        
    for shpath, shname in list(zip(local_paths_shell, shell_scripts)):
        
        ftp_client.put(localpath = shpath, remotepath = remote_path + shname, confirm = False)
        
        shell_remote_paths.append(remote_path + shname)
        
        
    ftp_client.close()
    
    # before submitting to queue system, we correct possible mistakes between DOS / UNIX
    for input_path in input_remote_paths:
        new_path = input_path[:-4]+'_new.inp'
        stdin, stdout, stderr = nodehost.exec_command('tr "\r" "\n" < '+input_path+' > '+new_path)
        stdin, stdout, stderr = nodehost.exec_command('mv '+new_path+' '+input_path)
        
    for shell_path in shell_remote_paths:
        new_path = shell_path[:-3]+'_new.sh'
        stdin, stdout, stderr = nodehost.exec_command('tr "\r" "\n" < '+shell_path+' > '+new_path)
        stdin, stdout, stderr = nodehost.exec_command('mv '+new_path+' '+shell_path)
        
    nodehost.close()
    host.close()
    
    return shell_remote_paths


def run_simulation_SGE(shell_remote_paths, server = 'mastercr1.igfae.usc.es', node = 'nodo014', \
                       username = 'sergio.cabana', manual = True):
    
    ''' Submits to SGE queue system the shell_scripts prepared with mazepin.setup_simulation_SGE()
    
        ==================================================================================
        
        shell_remote_paths (list of str): list of paths to shell scripts inside remote machine. Output of setup_simulation_SGE
        
        server, node, username (str): same as in setup_simulation_SGE
        
        manual (bool): If True, asks for confirmation to submit every script
    
    '''
    
    # we connect to our remote host
    print('PASSWORD INPUT: Be careful, it might be visible in some terminals')
    print('If you did not get a warning before this, it should work properly')
    
    password = getpass('Please introduce your remote machine password: ')
    
    host = paramiko.SSHClient()
    host.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    host.connect(server, username=username, password=password)

    # we open a communication channel
    hosttransport = host.get_transport()
    dest_addr = (node, 22)
    local_addr = (server, 22)
    hostchannel = hosttransport.open_channel("direct-tcpip", dest_addr, local_addr)

    # connection with node 
    nodehost = paramiko.SSHClient()
    nodehost.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    nodehost.connect(node, username=username, password=password, sock=hostchannel)
    
    # now we submit to queque system:
    i = 0
    for shell_path in shell_remote_paths:
        if manual:
            i += 1
            input(f'Press Enter to submit bash script {i}: ')
            print('Done')
        stdin, stdout, stderr = nodehost.exec_command('qsub ' + shell_path)
        
    # we check that everything is working
    
    print('Everything should be submitted by now. Let\'s see what qstat returns (Wait for 5 secs):')
    print('===========================================================\n')
    
    time.sleep(5)
    stdin, stdout, stderr = nodehost.exec_command('qstat')
    for line in stdout.readlines():
        print(line)
        
    nodehost.close()
    host.close()
    

# def simulator(task_names, basics, trajects, sim_controls, exports, extras, jobID, \
#               RASPASS = False, upgoing = False, ZHAireS = False, \
#               ZHAireS_control = [], remote_main = '/home2/', user = 'sergio.cabana/', remote_dir = '', \
#               exe = ['aires/bin/ZHAireSRASPASS', 'SpecialPrimaries/RASPASSprimary', 'SpecialPrimaries/uprimary'], \
#               program = 'ZHAireSRASPASS', autorename = True, local_savepath = '', \
#               eco = False, server = 'mastercr1.igfae.usc.es', node = 'nodo014', \
#               username = 'sergio.cabana'):
    
#     ''' Full simulation process.
    
#         Creates input files, shell scripts, submits to remote machine and runs simulations
        
#         ============================================================================
        
#         task_names (list): names of tasks (max 64 char), one per each input file required
        
#         basics (list): set of basic parameters (arg of mazepin.input_file()). one per each input file
        
#         trajects (list): set of trajectory params (arg of mazepin.input_file()). one per each input file
        
#         sim_control (list): set of parameters to control simulation (arg of mazepin.input_file()). Default options at mazepin_aux. one per each input file
        
#         exports (list): lists of [tab, part, option] to export in each input (arg of mazepin.input_file()). one per each input file.
        
#         extras (list): lists of extra instructions to include in each input file (arg of mazepin.input_file()). one per each input file.
        
#         jobID (str): Identifier of remote job (qsub)
        
#         RASPASS, upgoing (bool): Special particle options
        
#         ZHAireS (bool): input files include radio emission calculations
        
#         ZHAireS_control: (list): set of parameters to run ZHAireS. (arg of mazepin.input_file()).
        
#         remote_main (str, ends with /): root directory in remote machine
        
#         user (str, edns with /): remote machine user. All directories in remote machine should be like /remote_dir/user/something
        
#         remote_dir (str, ends with /): directory inside remote_dir/user/ to upload files and store outputs
            
#         exe (list of str): list of paths inside remote_dir/user/ where executables to run special particles are
        
#         program: name of executable that runs the simulations
        
#         autorename: bool to include loops in shell script to rename tables with opt a

#         local_savepath (str): directory to save created files
        
#         eco (bool): True to send only one shell script to remote (otherwise, one per input file to parallelize)
        
#         server: host remote machine
        
#         node : where to ssh to from remote server
        
#         username: name of remote user
#         ===========================================================================================
        
#         Though there are quite a few inputs, most of themo wont change between input files (usually,
#         only one or two parameters change, mostly energies or trajectories)
        
#         Will ask for remote machine password (assumes same password for both remote server and node)
#     '''
#     lengths = [len(task_names), len(basics), len(trajects), len(sim_controls), len(exports), len(extras)]
    
#     if not all([l == lengths[0] for l in lengths]):
#         raise TypeError('Mismatch of lengths of task_names, basics, sim_controls, exports, extras')
        
#     if ZHAireS and len(ZHAireS_control) != len(task_names):
#         raise TypeError('Length of ZHAireS_control does not match number of input files')
        
    
#     # We create all the input files
#     simulations = list(zip(task_names, basics, trajects, sim_controls, exports, extras))
#     counter     = 0
    
#     local_paths = []
#     input_files = []
    
#     for taskname, basic, traj, sim_control, exp, ext in simulations:
        
#         ZHScon = ZHAireS_control[counter] if ZHAireS else []
        
#         path, name = input_file(task_name = taskname, basic = basic, trajectory = traj, \
#                                 sim_control = sim_control, RASPASS = RASPASS, upgoing = upgoing, \
#                                 ZHAireS = ZHAireS, ZHAireS_control = ZHScon, exports = exp, \
#                                 extra = ext, save_path = local_savepath)
            
#         local_paths.append(path)
#         input_files.append(name)
        
#         counter += 1
        
#     # we create all the required shell_scripts
    
#     local_paths_shell = []
#     shell_scripts     = []
    
#     if eco:
#         path, name = shell_script(job_name = jobID, input_file_dir = remote_dir, \
#                                   input_files = input_files, local_dir = local_savepath, \
#                                   main = remote_main, user = user, exe = exe, program = program,\
#                                   autorename = autorename)
        
#         local_paths_shell.append(path)
#         shell_scripts.append(name)
            
#     else:
#         for i in range(len(input_files)):
#             path, name = shell_script(job_name = jobID+str(i), input_file_dir = remote_dir, \
#                                       input_files = [input_files[i]], local_dir = local_savepath, \
#                                       main = remote_main, user = user, exe = exe, program = program, \
#                                       autorename = autorename)
            
#             local_paths_shell.append(path)
#             shell_scripts.append(name) 
            
#     # we connect to our remote host
    
#     password = str(input('Please introduce your remote machine password (careful, it is visible): '))
    
#     host = paramiko.SSHClient()
#     host.set_missing_host_key_policy(paramiko.AutoAddPolicy())
#     host.connect(server, username=username, password=password)

#     # we open a communication channel
#     hosttransport = host.get_transport()
#     dest_addr = (node, 22)
#     local_addr = (server, 22)
#     hostchannel = hosttransport.open_channel("direct-tcpip", dest_addr, local_addr)

#     # connection with node 
#     nodehost = paramiko.SSHClient()
#     nodehost.set_missing_host_key_policy(paramiko.AutoAddPolicy())
#     nodehost.connect(node, username=username, password=password, sock=hostchannel)

#     # open sftp channel to upload files
    
#     ftp_client  = host.open_sftp()
    
#     remote_path = remote_main + user + remote_dir
    
#     input_remote_paths = []
    
#     shell_remote_paths = []
    
#     for inppath, inpname in list(zip(local_paths, input_files)):
        
#         ftp_client.put(localpath = inppath, remotepath = remote_path + inpname, confirm = False)
        
#         input_remote_paths.append(remote_path + inpname)
        
#     for shpath, shname in list(zip(local_paths_shell, shell_scripts)):
        
#         ftp_client.put(localpath = shpath, remotepath = remote_path + shname, confirm = False)
        
#         shell_remote_paths.append(remote_path + shname)
        
        
#     ftp_client.close()
    
#     # before submitting to queue system, we correct possible mistakes between DOS / UNIX
#     for input_path in input_remote_paths:
#         new_path = input_path[:-4]+'_new.inp'
#         stdin, stdout, stderr = nodehost.exec_command('tr "\r" "\n" < '+input_path+' > '+new_path)
#         stdin, stdout, stderr = nodehost.exec_command('mv '+new_path+' '+input_path)
        
#     for shell_path in shell_remote_paths:
#         new_path = shell_path[:-3]+'_new.sh'
#         stdin, stdout, stderr = nodehost.exec_command('tr "\r" "\n" < '+shell_path+' > '+new_path)
#         stdin, stdout, stderr = nodehost.exec_command('mv '+new_path+' '+shell_path)
        
#     # now we submit to queque system:

#     for shell_path in shell_remote_paths:
#         stdin, stdout, stderr = nodehost.exec_command('qsub ' + shell_path)
    
#     # we check that everything is working
    
#     print('Everything should be submitted by now. Let\'s see what qstat returns (Wait for 5 secs):')
#     print('===========================================================\n')
    
#     time.sleep(5)
#     stdin, stdout, stderr = nodehost.exec_command('qstat')
#     for line in stdout.readlines():
#         print(line)
        
        
#     nodehost.close()
#     host.close()

# def Xs_to_dist_old(X, RD, RH, theta, prec = .05):
#     ''' Converts slanted depth [g/cm2] to distance covered along shower axis,
#         for a RASPASS trajectory. VERY approximate result 
#         (rough and easy numerical integration)
    
#         X: values of slanted depth [g/cm2] (numpy.ndarray)
        
#         RD: RASPASSDistance [km]
    
#         RH: RASPASSHeight [km]
    
#         theta: PrimaryZenAngle [deg]
    
#         prec: precission for the result. Default 50m
#     '''

#     thetarad = theta *np.pi/180.
    
#     def integrand(h):
#         return rho_Lin(h)/np.sqrt(1-((RT+RH)/(RT+h)*np.sin(thetarad))**2)
#         # rho_Lin is the Linsley atmosphere, defined in mazepin_aux
        
#     X      = np.sort(X)    # we sort values (just in case)
#     X_last = max(X)        # last value, will be the highest
    
#     hmin   = (RT+RH)*np.sin(thetarad) - RT
#     hmin   = hmin if hmin > 0 else 0  #minimum height in atmosphere
    
#     x      = 0.
#     L      = 0.
    
#     dist   = []
#     index  = 0
    
#     while x < X_last:
#         h1    = h_RAS(L, RD, RH, theta)
#         h2    = h_RAS(L+prec, RD, RH, theta)
        
#         if hmin < h1 or hmin > h2:
#             delta_x = abs(integrand(.5*(h1+h2))*(h2-h1))*1e5
#         else:
#             delta_x = abs(integrand(.5*(hmin+h1))*(hmin-h1))*1e5 + \
#                       abs(integrand(.5*(h2+hmin))*(h2-hmin))*1e5
        
#         # we find matter traversed when we move a distance prec
#         # absolute value to deal with increasing/decreasing heights
#         # take into account interval where height starts increasing
        
#         x_new = x + delta_x

#         if x_new > X[index]: 
#             # X[index] is between x and x_new = x+delta_x
#             # we use a linear interpolation between the last two values
            
#             distance = L + prec/delta_x * (X[index]-x)
#             index   += 1
#             dist.append(distance)
        
#         x  = x_new
#         L += prec
        
#     return np.array(dist)

