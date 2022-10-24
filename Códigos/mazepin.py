############################## MAZEPIN v0.7.1 #################################
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
from   scipy.interpolate import interp1d
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
    
        RD: RASPASSDistance [km]
        
        RH: RASSPASSHeight  [km]
        
        theta: PrimaryZenAngle (deg)
    '''
        
    thetarad = theta *np.pi/180.
    
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
    
    v = h_IP(RD, RH, theta)
    c = cos_localtheta(v, theta, RASPASSHeight = RH)
    
    return np.abs(np.sqrt((RT+v)**2+L**2-2*L*(RT+v)*c)-RT)

def dist_to_Xs(dist, RD, RH, theta, step = .025):
    ''' Converts distance [km] to slanted depth [g/cm2] covered along shower 
        axis, for a RASPASS trajectory. VERY approximate result 
        (rough and easy numerical integration)
    
        d: values of distance [km] (numpy.ndarray)
        
        RD: RASPASSDistance [km]
    
        RH: RASPASSHeight [km]
    
        theta: PrimaryZenAngle [deg]
        
        step: Step of integration. Default 25 m
        
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
       

def atmos_size(RH, theta, atmos_height = 113, stop = 'atmos'):
    ''' Returns total distance and traversed matter for a RASPASS trajectory
        starting at the top ot the atmosphere (fixes RASPASSDistance so that 
        happens)
        
        RH : RASPASSHeight [km]
        
        theta : PrimaryZenAngle [deg]
        
        atmos_height : height of the atmosphere [km] and starting
            point of the trajectory (default 113 km)
            
        stop : point up to which we measure distance and traversed matter
            'atmos' : top of the atmosphere
            'zaxis' : crossing of shower axis with the vertical of observer
            'hmin'  : minimum height reached in the atmosphere
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
        
        
def Xs_to_dist(X, RD, RH, theta, atmos_height = 113):
    ''' Converts slanted depth [g/cm2] to distance covered along shower axis,
        for a RASPASS trajectory. VERY approximate result 
        (rough and easy numerical integration). Inverts previous function.
    
        X: values of slanted depth [g/cm2] (numpy.ndarray)
        
        RD: RASPASSDistance [km]
    
        RH: RASPASSHeight [km]
    
        theta: PrimaryZenAngle [deg]

    '''
    
    RD_top = RAS_Dist(atmos_height, theta, RH = RH)
    
    offset_dist = RD - RD_top if RD - RD_top > 0 else 0 # distance outside atmosphere

    L_full, X_full = atmos_size(RH, theta)
    
    d = np.linspace(0, L_full+offset_dist, 1000)
    
    Xs = dist_to_Xs(d, RD, RH, theta)
    
    inverse = interp1d(Xs, d)
    
    return np.array([inverse(x) if x < X_full else L_full+offset_dist for x in X])
    
def codes():
    return dict_tab, dict_part

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

def pathfinder(rootdir, tabs, part, sim = [''], sep = [''], verbose = False):
    ''' Returns the list of paths inside rootdir corresponding to the requested
        tables.
        
        tabs (list of int): types of tables. See mazepin.table_finder() .
        
        part (list of int): types of particles. See mazepin.table_finder() .
        
        sim (list of str): constraints for file names (all selected files must 
        these strings). Default '', all files with the right extension are kept.
        I am using my rules for naming files (see notes at GitHub repo). 
        
        sep (list of str): Special distinctions. Output will indicate which files
        contain these special strings. Default '', all files with the right 
        extension are kept
        
        Output (with separators sep = [sep1, sep2]):
            
        [[tab1, part1, sep1, path1], [tab1, part1, sep2, path2], 
         [tab1, part2, sep1, path3], [tab1, part2, sep2, path4],
         [tab2, part1, sep1, path5], ...]
        
        ** sorry for the spaghetti, I know it is a huge mess :_( **
    '''
    paths = []
    
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:     # loop over files in rootdir
            for tab in tabs:   # loop over requested tables
                for p in part: # loop over requested particles
                
                    if tab != 16:  # NORMAL AIRES TABLES (JUST LOOP OVER FILES AND CHECK EXTENSIONS)
                    
                        table_ext = '.t'+str(tables[p,tab]) 
                        # extension of table tab, particle p (table at mazepin_aux)
                    
                        if file.endswith(table_ext) and all([c in file for c in sim]):
                        # if our file has the right extension and all constraints
                            for s in sep: # loop over separations
                                if s == '' or s in file: #we take into account distinctions
                                    paths.append([str(tab), str(p), s, subdir + os.sep + file])
                                    # save type of table, particle, separation and path
                                    
                    elif tab == 16: # energy per particle, we need 2 tables to calculate
                    
                        table_ext1 = '.t'+str(tables[p,0]) 
                        table_ext2 = '.t'+str(tables[p,2])
                        # extension of tables 0 and 2, particle p (table at mazepin_aux)
                        
                        if file.endswith(table_ext1) and all([c in file for c in sim]):
                        # if our file has the right extension and all constraints
                            for s in sep: # loop over separations
                                if s == '' or s in file: #we take into account distinctions
                                    filename = file[:-6] 
                                    #name of task, should have both extensions inside rootdir
                                    paths.append([str(tab), str(p), s, \
                                                 [subdir + os.sep + file, subdir + os.sep + filename+table_ext2]])
                                    # save both paths in the same element
                                    
    if len(sep) > 1 and len(paths) != len(part) * len(sep) * len(tabs):
        raise TypeError('Some tables are not available inside the main directory')
        
    order = []
    for p in part:
        for s in sep:
            for elem in paths:
                if str(p) in elem and s in elem:
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
               RASPASS = False, Distance = False, first_int = False, \
               traject = []):
    ''' Returns data to be plotted, and a collection of suitable labels
    
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
        
        if error_type == 'sigma':
            err = data[:,3]
        elif error_type == 'RMS':
            err = data[:,2]
        else:
            raise TypeError('Please introduce a valid error_type ("sigma" or "RMS"")')
        
    else:
        data1, grd = readfile(file[0]) # longitudinal development
        data2, _  = readfile(file[1])  # energy longitudinal development data
        
        xdata = data1[:,0]; y1 = data1[:, 1]; y2 = data2[:,1]
        ydata = [a/b if b!=0 else 0 for a, b in list(zip(y2, y1))] # divide both datasets
        
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
            raise TypeError('Please introduce a valid error_type ("sigma" or "RMS"")')
        
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
        
        RD_top       = RAS_Dist(113, theta, RH = 0) # RASPASSDistance of an inj ponit at atmosphere limit
        traj_input   = True
        
        
    ######################### X AXIS CONVERSIONS ##########################
    if not RASPASS:  # "normal" upgoing or downgoing showers, RH = 0
        
        if Distance and tab in tabs_x_depth: # we want to convert x_data to distances along axis 
        
            x_axis_label = r'Dist. along axis [km]' 
            
            if not traj_input:
                raise TypeError('Trajectory parameters are needed')
                
            if not slant:
                raise TypeError('Conversor Xs_to_dist works with X_slanted as input')

            grd_km = Xs_to_dist(np.array([grd]), RD_top, 0, theta) if grd != 'None' else RD_top
            # position of ground from atmosphere limit, along axis
            
            xdata  = Xs_to_dist(xdata, RD_top, 0, theta)
            
            if first_int: # we start counting at the injection height, upwards or downwards
                xdata = np.array([grd_km - x - RD for x in xdata]) if UG else np.array([x - RD for x in xdata])
            else:
                xdata = np.array([grd_km - x for x in xdata]) if UG else np.array(xdata)
                
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
                raise TypeError('RASPASS tables by default export distance in x axis (without any option). Do not use slant')
            
        elif slant and not Distance and tab in tabs_x_depth:
            
            x_axis_label = r'$X_s$ [$\mathrm{g/cm^2}$]' 
            

    return [xdata, ydata, err, label], x_axis_label, y_axis_label


def Aires_Plot(input_data, error_type = 'sigma', UG = False, slant = False, \
               RASPASS = False, Distance = False, first_int = False, \
               trajects = [], xlim = [], ylim = [], xscale = 'linear', \
               yscale = 'linear', autolabel = True, graph_type = 'step', labels = [], \
               size = 12, legend = True, title = '', loc_leg = 'best', \
               fmt = '-', marker = 'o', linewidth = 2.):
    
    
    ''' Plots AIRES outputs, preprocessed using pathfinder and Aires_data
    
        input_data : output of pathfinder. The selection of graphs is done there
        
        Allows to plot different tables in the same canvas, as long as it is possible
        (i.e., same X_axis coordinate)
        
        error_type, UG, slant, RASPASS, Distance, first_int, trajects
        are Aires_data options
        
        If traject is not needed, place an empty list [] for each plot
        
        xlim, ylim (list): limits for axes
        
        xscale, yscale ("linear", "log"): scales for axes
        
        autolabel (bool): Sets automatic labels for each plot
        
        graph_type: Valid types are "step" or "errorbar"
        
        labels: if autolabel is False, list of labels by hand (one for plot,
                in the order of input_data)
        
        size : size of text
        
        legend (bool): sets legend in graph
        
        title (str): title of graph
        
        loc_leg : position of legend
        
        fmt : if graph_style == errorbar, format of plot
        
        marker : if graph_style == errorbar, marker of points
        
        linewidth : for both types of graphs
        
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
    
    fig = plt.figure()

    ax  = fig.add_subplot(111)
    plt.xticks(fontsize = size)
    plt.yticks(fontsize = size)
    
    # scientific notation for offset text

    form = mpl.ticker.ScalarFormatter(useMathText=True) 
    form.set_scientific(True) 
    form.set_powerlimits((-1,1))
    ax.yaxis.set_major_formatter(form)
    
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
        
        if autolabel:
            if multi_y:
                glabel = dict_tab_yleg[tabs[counter]]+', '+lab # dict in mazepin_aux
            else:
                glabel = lab
        else:
            glabel = labels[counter]
            
        if graph_type == 'step':
            ax.step(x, y, linewidth = linewidth, where = 'mid', label = glabel)

        elif graph_type == 'errorbar':
            ax.errorbar(x, y, err, fmt = fmt, marker=marker, linewidth=linewidth, label = glabel)

        else:
            raise TypeError('Valid graph styles are "step" and "errorbar"')
            
        counter += 1
    
    ax.set_xlabel(x_axis_label, size = size)
    ax.set_ylabel(y_axis_label, size = size)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.set_title(title, size = size)
    
    text = ax.yaxis.get_offset_text()
    text.set_size(size)
    
    if UG and not Distance and not slant:
        ax.invert_xaxis()
        
    if RASPASS and Distance:
        ax.invert_xaxis()
    
    if legend:
        ax.legend(loc = loc_leg, prop = {'size':size})
        
    if len(xlim) == 2:
        ax.set_xlim(xlim[0], xlim[1])
    
    if len(ylim) == 2:
        ax.set_ylim(ylim[0], ylim[1])
        
    return fig, ax
        

def traject_finder(files, RASPASS = False, UG = False):
    ''' Returns trajecory parameters given the name of the file
        **WARNING** It uses my personal conventions for naming files
        
        files: output of pathfinder
        
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
            inj_h = chars[2][:-2]
            theta = chars[3][:-3]
            
            trajects.append([float(inj_h), float(theta)])
    
    return trajects
            

def shell_script(job_name, input_file_dir, input_files, local_dir = '', \
                 main = '/home2/', user = 'sergio.cabana/', \
                 exe = ['aires/bin/ZHAireSRASPASS', 'SpecialPrimaries/RASPASSprimary'],\
                 program = 'ZHAireSRASPASS'):
    
    ''' Creates a shell script (.sh extension) that can be submitted to queue 
        systems (at least inside IGFAE nodes, with scratch). 
        It is adapted for the cluster I use but easy to modify
        
        job_name: job ID for queue system
        
        local_dir: directory to save the created .sh file
        
        input_file_dir: directory of input files (inside main/user/). Ends with /
        
        input_files : list of names of input files
        
        main : default root directory
        
        user : directory of current user
        
        exe  : paths (inside main) to executable binaries to import and use

        program : name of binary executable to run the input files
    '''

    path = os.path.join(local_dir, 'script_shell_'+job_name+'.sh')
    
    header = (        
r'''############ Created with mazepin.shell_script() ################

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
    
    scratch_dir = r'/scratch/'+user+'Aires_tmp/'+job_name+'/'

    with open(path, 'w') as f:
        
        f.write(header)
        f.write(r'[[ -d '+ scratch_dir + ' ]] || mkdir -p '+ scratch_dir+'\n')
        f.write(r'cd ' + scratch_dir+'\n')
        
        f.write('############ EXECUTABLES ###########\n')
        for e in exe:
            
            f.write(r'cp ' + main + user + e + ' ' + scratch_dir+'\n')
        
        
        f.write('############ INPUT FILES ###########\n')
        for inp in input_files:
            
            f.write(r'cp ' + main + user + input_file_dir + inp + ' ' + scratch_dir+'\n')
            
        
        f.write('############ EXECUTION ###########\n')
        for inp in input_files:
            
            f.write(r'./' + program + ' < ' + inp + ' > ' + inp[:-4] + '.out'+'\n')
            
        
        f.write('############ OUTPUT SAVING AND EXIT ########### \n')

        f.write(r'cp *.* '+ main + user + input_file_dir+'\n')
        f.write(r'rm -rf ' + scratch_dir+'\n')
                    
        f.close()


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
