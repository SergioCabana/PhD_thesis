''' MAZEPIN_AUX
    Here I will collect some definitions of arrays, lists and dictionaries 
    that will be imported and used by the main module
'''
import numpy as np
import matplotlib as mpl

# colorblind friendly sequence of colors, at least for me ;)

colors = ['k', 'royalblue', 'r', 'gold', 'limegreen', 'navy', 'crimson', \
           'turquoise', 'darkorange', 'darkgreen']
    

# Linsley model for the atmosphere

b1, c1 = 1222.6562, 9941.8638e2
b2, c2 = 1144.9069, 8781.5355e2
b3, c3 = 1305.5948, 6361.4304e2
b4, c4 = 540.1778, 7721.7016e2
b5, c5 = 1, 1e9                  

# coefficients b [g/cm2], c [cm]
# see AIRES manual, table 2.1

def rho_Lin(h):
    if 0 <= h < 4:
        return b1*np.exp(-h*1e5/c1)/c1
    elif 4 <= h < 10:
        return b2*np.exp(-h*1e5/c2)/c2
    elif 10 <= h < 40:
        return b3*np.exp(-h*1e5/c3)/c3
    elif 40 <= h < 100:
        return b4*np.exp(-h*1e5/c4)/c4
    elif 100 <= h < 113:
        return b5/c5
    else:
        return 0
        
    
    
# dictionary of tables in AIRES 19.04.08 (table 16 is added by me)

dict_tab = {'0' : 'Longitudinal development',
            '1' : 'Unweighted longitudinal development',
            '2' : 'Energy longitudinal development',
            '3' : 'Lateral distribution',
            '4' : 'Unweighted lateral distribution',
            '5' : 'Energy distribution at ground',
            '6' : 'Unweighted energy distribution',
            '7' : 'Mean arrival time distribution',
            '8' : 'Number and energy of particles at ground vs shower number',
            '9' : 'Number of created particles',
            '10': 'Number of created entries',
            '11': 'Energy of created particles',
            '12': 'Longitudinal development of low energy particles',
            '13': 'Unweighted longitudinal development of low energy particles',
            '14': 'Energy longitudinal development of low energy particles',
            '15': 'Longitudinal development of deposited energy',
            '16': 'Energy (per particle) longitudinal development'}

# set of tables where x_axis is a traversed depth:
    
tabs_x_depth = ['0', '1', '2', '9', '10', '11', '12', '13', '14', '15', '16']

# set of tables where x_axis is energy:
    
tabs_x_ener = ['5', '6']

# set of tables where x_axis is distance to core:
    
tabs_x_core = ['3', '4', '7']

# dictionary of x_labels for plots of each table (default export specifications)

dict_tab_xlab = {'0' : r'$X_v$ [$\mathrm{g/cm^2}$]',
                 '1' : r'$X_v$ [$\mathrm{g/cm^2}$]',
                 '2' : r'$X_v$ [$\mathrm{g/cm^2}$]',
                 '3' : r'Distance to core [$\mathrm{m}$]',
                 '4' : r'Distance to core [$\mathrm{m}$]',
                 '5' : r'$E$ [$\mathrm{GeV}$]',
                 '6' : r'$E$ [$\mathrm{GeV}$]',
                 '7' : r'Distance to core [$\mathrm{m}$]',
                 '8' : 'Not implemented',
                 '9' : r'$X_v$ [$\mathrm{g/cm^2}$]',
                 '10': r'$X_v$ [$\mathrm{g/cm^2}$]',
                 '11': r'$X_v$ [$\mathrm{g/cm^2}$]',
                 '12': r'$X_v$ [$\mathrm{g/cm^2}$]',
                 '13': r'$X_v$ [$\mathrm{g/cm^2}$]',
                 '14': r'$X_v$ [$\mathrm{g/cm^2}$]',
                 '15': r'$X_v$ [$\mathrm{g/cm^2}$]',
                 '16': r'$X_v$ [$\mathrm{g/cm^2}$]'}

# dictionary of y_labels for plots of each table (default export specifications)

dict_tab_ylab = {'0' : r'Number of particles',
                 '1' : r'Number of particles',
                 '2' : r'$E$ [$\mathrm{GeV}$]',
                 '3' : r'Number of particles',
                 '4' : r'Number of particles',
                 '5' : r'Number of particles',
                 '6' : r'Number of particles',
                 '7' : r'Arrival time [$\mathrm{ns}$]',
                 '8' : 'Not implemented',
                 '9' : r'Number of created particles',
                 '10': r'Number of created entries',
                 '11': r'$E$ of created particles [$\mathrm{GeV}$]',
                 '12': r'Number of low $E$ particles',
                 '13': r'Number of low $E$ particles',
                 '14': r'$E$ of low $E$ particles [$\mathrm{GeV}$]' ,
                 '15': r'Deposited energy [$\mathrm{GeV}$]',
                 '16': r'$E$ (per particle) [$\mathrm{GeV}$]'}

# dictionary of y_labels for legends

dict_tab_yleg = {'0' : r'#',
                 '1' : r'#',
                 '2' : r'$E$',
                 '3' : r'#',
                 '4' : r'#',
                 '5' : r'#',
                 '6' : r'#',
                 '7' : r'Arrival time',
                 '8' : 'Not implemented',
                 '9' : r'# (created)',
                 '10': r'# (entries)',
                 '11': r'$E$ (created)',
                 '12': r'# (low E)',
                 '13': r'# (low E)',
                 '14': r'$E$ (low $E$)' ,
                 '15': r'Dep. $E$',
                 '16': r'$E$ (per part.)'}

# dictionary of particles in AIRES 19.04.08

dict_part = {'0' : 'gammas',
             '1' : 'e-',
             '2' : 'e+',
             '3' : 'mu+',
             '4' : 'mu-',
             '5' : 'pi+',
             '6' : 'pi-',
             '7' : 'K+',
             '8' : 'K-',
             '9' : 'n',
             '10': 'p',
             '11': 'anti-p',
             '12': 'Nuclei',
             '13': 'Other charged pcles.',
             '14': 'Other neutral pcles.',
             '15': 'e+e-',
             '16': 'mu+mu-',
             '17': 'pi+pi-',
             '18': 'K+K-',
             '19': 'All charged pcles.',
             '20': 'All neutral pcles.',
             '21': 'All pcles.',
             '22': 'All neutrinos'}

# dictionary of particles in AIRES 19.04.08 (LaTeX)

dict_part_tex = {'0' : r'$\gamma$',
                 '1' : r'$e^-$',
                 '2' : r'$e^+$',
                 '3' : r'$\mu^+$',
                 '4' : r'$\mu^-$',
                 '5' : r'$\pi^+$',
                 '6' : r'$\pi^-$',
                 '7' : r'$K^+$',
                 '8' : r'$K^-$',
                 '9' : r'$n$',
                 '10': r'$p$',
                 '11': r'$\bar{p}$',
                 '12': r'Nuclei',
                 '13': r'Other charged pcles.',
                 '14': r'Other neutral pcles.',
                 '15': r'$e^+e^-$',
                 '16': r'$\mu^+\mu^-$',
                 '17': r'$\pi^+\pi^-$',
                 '18': r'$K^+K^-$',
                 '19': r'All charged pcles.',
                 '20': r'All neutral pcles.',
                 '21': r'All pcles.',
                 '22': r'All neutrinos'}

# matrix of tables (see table_finder). You can generate this array with the
# script 'Table_index_creator', uploaded at GitHub

tables = np.array([[1001, 1301, 1501, 2001, 2301, 2501, 2801, 3001, 5001, 6001,\
                    6301, 6501, 7001, 7301, 7501, 7801],
       [1005, 1305, 1505, 2005, 2305, 2505, 2805, 9999, 5005, 6005, 6305,\
        6505, 7005, 7305, 7505, 7805],
       [1006, 1306, 1506, 2006, 2306, 2506, 2806, 9999, 5006, 6006, 6306,\
        6506, 7006, 7306, 7506, 7806],
       [1007, 1307, 1507, 2007, 2307, 2507, 2807, 9999, 5007, 6007, 6307,\
        6507, 7007, 7307, 7507, 7807],
       [1008, 1308, 1508, 2008, 2308, 2508, 2808, 9999, 5008, 6008, 6308,\
        6508, 7008, 7308, 7508, 7808],
       [1011, 1311, 1511, 2011, 2311, 2511, 2811, 9999, 5011, 6011, 6311,\
        6511, 9999, 9999, 9999, 9999],
       [1012, 1312, 1512, 2012, 2312, 2512, 2812, 9999, 5012, 6012, 6312,\
        6512, 9999, 9999, 9999, 9999],
       [1013, 1313, 1513, 2013, 2313, 2513, 2813, 9999, 5013, 6013, 6313,\
        6513, 9999, 9999, 9999, 9999],
       [1014, 1314, 1514, 2014, 2314, 2514, 2814, 9999, 5014, 6014, 6314,\
        6514, 9999, 9999, 9999, 9999],
       [1021, 1321, 1521, 2021, 2321, 2521, 2821, 9999, 5021, 6021, 6321,\
        6521, 9999, 9999, 9999, 9999],
       [1022, 1322, 1522, 2022, 2322, 2522, 2822, 9999, 5022, 6022, 6322,\
        6522, 9999, 9999, 9999, 9999],
       [1023, 1323, 1523, 2023, 2323, 2523, 2823, 9999, 5023, 6023, 6323,\
        6523, 9999, 9999, 9999, 9999],
       [1041, 1341, 1541, 2041, 2341, 2541, 2841, 9999, 5041, 6041, 6341,\
        6541, 9999, 9999, 9999, 9999],
       [1091, 1391, 1591, 2091, 2391, 2591, 2891, 9999, 5091, 6091, 6391,\
        6591, 7091, 7391, 7591, 7891],
       [1092, 1392, 1592, 2092, 2392, 2592, 2892, 3092, 5092, 6092, 6392,\
        6592, 7092, 7392, 7592, 7892],
       [1205, 1405, 1705, 2205, 2405, 2705, 2905, 3205, 5205, 6205, 6405,\
        6705, 7205, 7405, 7705, 7905],
       [1207, 1407, 1707, 2207, 2407, 2707, 2907, 3207, 5207, 6207, 6407,\
        6707, 7207, 7407, 7707, 7907],
       [1211, 1411, 1711, 2211, 2411, 2711, 2911, 9999, 5211, 6211, 6411,\
        6711, 9999, 9999, 9999, 9999],
       [1213, 1413, 1713, 2213, 2413, 2713, 2913, 9999, 5213, 6213, 6413,\
        6713, 9999, 9999, 9999, 9999],
       [1291, 1491, 1791, 2291, 2491, 2791, 2991, 3291, 5291, 6291, 6491,\
        6791, 7291, 7491, 7791, 7991],
       [1292, 1492, 1792, 2292, 2492, 2792, 2992, 3292, 5292, 6292, 6492,\
        6792, 7292, 7492, 7792, 7992],
       [1293, 1493, 1793, 2293, 2493, 2793, 2993, 3293, 5293, 6293, 6493,\
        6793, 7293, 7493, 7793, 7993],
       [9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 6296, 6496,\
        6796, 9999, 9999, 9999, 9999]])
    

dict_basic_IDL = {'0': 'PrimaryParticle',
                  '1': 'PrimaryEnergy',
                  '2': 'GeomagneticField',
                  '3': 'PropagatePrimary',
                  '4': 'Site',
                  '5': 'GroundAltitude'}

dict_control_IDL = {'0': 'TotalShowers',
                    '1': 'RunsPerProcess',
                    '2': 'ShowersPerRun',
                    '3': 'RandomSeed',
                    '4': 'ObservingLevels',
                    '5': 'ThinningEnergy',
                    '6': 'ThinningWFactor',
                    '7': 'SaveInFile grdpcles',
                    '8': 'SaveInFile lgtpcles',
                    '9': 'ElectronCutEnergy',
                    '10': 'ElectronRoughCut',
                    '11': 'GammaCutEnergy',
                    '12': 'GammaRoughCut'}

    
