import numpy as np
import matplotlib.pyplot as plt
import mazepin as maz

#%% mazepin.Aires_plot example

# I have uploaded a set of AIRES output tables to show how this works
# There I have included tables for normal, upgoing and RASPASS showers,
# exported with the default option (files ending in _vert) and with the 
# "Opt a" option (files ending with _slan).

# The downgoing and upgoing shower files were obtained with AIRES 2.8.4a,
# I truly hope nothing has changed in the matter of output format :_| 

rootdir = 'mazepin_examples_files'

#%%    ################ NORMAL (DOWNGOING) SHOWERS #########################

# Uploaded files: long. dev. of e+e-, varying theta

# remember the codes used:
table_codes, particle_codes = maz.codes()
print(table_codes)
print('\n=======================\n')
print(particle_codes)

# ask for longitudinal development
tabs = [0]
# ask for e+e-
part = [15]
# constraint the simulation files: ask for strings that all file names must have
sim = ['dg', 'prot', '1e18eV', '100km', '00pdeg', 'vert']
# separations: give string to plot separately the datasets corresponding to each string
sep = ['00deg', '50deg', '85deg']

# retrieve the files
files = maz.pathfinder(rootdir, tabs, part, sim, sep)
# get the trajectories
trajects = maz.traject_finder(files)
# and plot
fig, ax = maz.Aires_Plot(files, trajects = trajects, graph_type = 'step')

# Let's ask for slanted data. We just have to modify sim
sim = ['dg', 'prot', '1e18eV', '100km', '00pdeg', 'slan']

# retrieve the files again
files = maz.pathfinder(rootdir, tabs, part, sim, sep)
# get the trajectories
trajects = maz.traject_finder(files)
# and plot, using slant = True
fig2, ax2 = maz.Aires_Plot(files, slant = True, trajects = trajects, graph_type = 'step')

# Let's now ask for distance in the x axis. Conversor works with slanted data,
# so our current files will do the work. Keep slant = True
fig3, ax3 = maz.Aires_Plot(files, slant = True, Distance = True, trajects = trajects, graph_type = 'step')

# I you have data averaged over several showers, you can also plot errorbars 
fig4, ax4 = maz.Aires_Plot(files, error_type = 'sigma', \
                           slant = True, Distance = True, trajects = trajects, graph_type = 'errorbar')

#%% ######################### UPGOING SHOWERS ###########################
# Uploaded files: long. dev. of e+e-, varying theta and injection height

# we start the same way

tabs = [0]
part = [15]
sim = ['upr', 'prot', '1e18eV', '05km', '00pdeg', 'vert'] # Injection height at 5km !
sep = ['95deg', '100deg', '110deg', '120deg'] # we separate different zenith angles

# retrieve the files
files = maz.pathfinder(rootdir, tabs, part, sim, sep)
# get the trajectories
trajects = maz.traject_finder(files, UG = True) # use UG = True to send angles to 0-90 deg

# and plot, using UG = True (it turns the x axis, so you start in higher Xv and go upper in the atmosphere)
fig, ax = maz.Aires_Plot(files, trajects = trajects, graph_type = 'step', UG = True)


# Let's ask for slanted data. We just have to modify sim
sim = ['upr', 'prot', '1e18eV', '05km', '00pdeg', 'slan']

# retrieve the files again
files = maz.pathfinder(rootdir, tabs, part, sim, sep)
# get the trajectories
trajects = maz.traject_finder(files, UG = True)
# and plot, using slant = True
fig2, ax2 = maz.Aires_Plot(files, slant = True, UG = True, trajects = trajects, graph_type = 'step')

# What if we do not start counting matter at the point of first interaction?
# It seems that e+e- reach maximum after different amounts of matter!
# Just a problem of offsets, be careful with that!

fig3, ax3 = maz.Aires_Plot(files, slant = True, UG = True, first_int = False,\
                           trajects = trajects, graph_type = 'step')
    
# Let's now ask for distance in the x axis. Conversor works with slanted data,
# so our current files will do the work. Keep slant = True
fig4, ax4 = maz.Aires_Plot(files, slant = True, Distance = True, UG = True,\
                            trajects = trajects, graph_type = 'step')

# If we do not count from the first interaction, something similar would happen
fig5, ax5 = maz.Aires_Plot(files, slant = True, Distance = True, UG = True, first_int = False, \
                            trajects = trajects, graph_type = 'step')

#%% ###################### RASPASS SHOWERS ##############################
# Uploaded files: Longitudinal development of number and energy of e+e-, mu+mu-;
# varying RASPASSHeight. RASPASSDistance was set to inject at 100 km of height

# we do the same, ask for table and particle, select files and separate them

# SUPER WARNING: When you do a RASPASS simulation, the default export option
# returns tables with DISTANCE on the x axis, and the "Opt a" option returns
# tables with slanted depth in the x axis. There is no interest in using vertical 
# traversed matter in these kind of geometry

# We do basically the same, and ask for slanted data
tabs = [0]
part = [15]
sim = ['RAS', 'gamma', '1e18eV', 'slan']
sep = ['15km', '20km', '25km', '30km', '35km', '40km']

#retrieve the files
files = maz.pathfinder(rootdir, tabs, part, sim, sep)

# get the trajectories, for RASPASS
trajects = maz.traject_finder(files, RASPASS = True)

# and plot
fig, ax = maz.Aires_Plot(files, RASPASS = True, slant = True, trajects = trajects, graph_type = 'step')

# now, we retrieve "_vert" files, exported with the default options and thus containing
# distances on the x axis

sim = ['RAS', 'gamma', '1e18eV', 'vert']
files = maz.pathfinder(rootdir, tabs, part, sim, sep)

# the other things (tab, part, sim, trajects) are the same

# plot using Distance = True
fig2, ax2 = maz.Aires_Plot(files, RASPASS = True, Distance = True, slant = True, \
                          trajects = trajects, graph_type = 'step')
    
# WRONG, DO NOT USE SLANT DATA (it's just for remembering, if you were to set slant = False
# and use slant data, the program would not notice so be awake). Comment the previous line please ;)

fig2, ax2 = maz.Aires_Plot(files, RASPASS = True, Distance = True, slant = False, \
                         trajects = trajects, graph_type = 'step')
    
    
# You can also set more than one particle and table per graph
tabs = [0]
part = [15, 16] # two particles
sim = ['RAS', 'gamma', '1e18eV', 'vert']
sep = ['20km', '30km', '40km']

files = maz.pathfinder(rootdir, tabs, part, sim, sep)
trajects = maz.traject_finder(files, RASPASS = True)

fig3, ax3 = maz.Aires_Plot(files, RASPASS = True, Distance = True, slant = False, \
                         trajects = trajects, graph_type = 'step', yscale = 'log')
    
tabs = [0, 2] # two tables
part = [16]
sim = ['RAS', 'gamma', '1e18eV', 'vert']
sep = ['20km', '30km', '40km']

files = maz.pathfinder(rootdir, tabs, part, sim, sep)
trajects = maz.traject_finder(files, RASPASS = True)

fig4, ax4 = maz.Aires_Plot(files, RASPASS = True, Distance = True, slant = False, \
                         trajects = trajects, graph_type = 'step')

#%% mazepin.simulator example

## This is commented to avoid problems, like running the whole script
 
# though there are a lot of inputs, I hope copy/paste will speed up all the preparation
# mazepin.simulator tries to avoid working with remote machines, creating and
# submitting directly all the necessary files to run the desired simulations

###################### task names for each input file ##########################

# task_names = ['input_file1', 
#               'input_file2']

# #################### basic parameters for each input file #######################

# # primary, energy in eV, GeomagneticField, PropagatePrimary, Site, Ground (with units, g/cm2 or m)

# # empty string means input file will not contain the IDL, and thus go with default

# basics = [
#           ['Proton', '1e15', 'On', 'Off', 'SouthPole', ''],
#           ['Proton', '1e15', 'On', 'Off', 'SouthPole', '']
#           ]

# #################### trajectory parameters############# 
# #In a RASPASS shower we have 5 params, for the rest 3 params

# trajects = [
#             ['100', '0', '0'],
#             ['100', '65', '0']
#             ] 
# # vertical downgoing and inclined downgoing

# ####################### simulation control basic parameters ##################

# #   TotalShowers, RunsPerProcess, ShowersPerRun, RandomSeed, ObservingLevels, ThinningEnergy,
# #   ThinningWFactor, SaveInFile grdpcles, SaveInFIle lgtpcles, ElectronCutEnergy,
# #   ElectronRoughCut, GammaCutEnergy, GammaRoughCut
# # Default options at mazepin_aux (default_sim_control and default_sim_control_ZHS for ZHAireS simulations)
 
# sim_controls = [
#                default_sim_control_ZHS,
#                ['1', 'Infinite', '1', '', '350', '1e-4 Relative', '0.06', 'None', 'None', '', '', '', '']
#                ]

# ######################## tables to export ###########################

# # in format [tab, part, option] with mazepin codes
# # Here, for the first input file, we export long. dev. of e+e-, normal and slanted
# # and for the second input file, the same and also long. dev of mu+mu-, and energy of mu+mu- (slanted)

# exports = [ 
#           [[0, 15, ''], [0, 15, 'a']],
#           [[0, 15, ''], [0, 15, 'a'], [0, 16, ''], [2, 16, 'a']] 
#           ]

# ###################### extra IDL instructions ##########################
# # not considered (not basic) for each input file

# extras = [
#          ['Atmosphere Linsley', 'PhotoNuclear Off'],
#          ['Atmosphere Linsley', 'PhotoNuclear Off']
#          ]

# ########################## ZHAireS control ###############################

# freqs   = [100, 300, 500]

# antenas = [
#           [0, 0, 0],
#           [0, 10, 0],
#           [10, 0, 0],
#           [0, -10, 0],
#           [-10, 0, 0],
#           ]

# direct_ZHS_inst = [
#                   ['RemoveZeroOutput On'],
#                   ['RemoveZeroOutput On']
#                   ]

# # On/Off for FresnelTime, On/Off for FresenlFreq, antena coordinates, frequency list, other direct instructions

# ZHAireS_control = [
#                   ['On', 'On', antenas, freqs, direct_ZHS_inst],
#                   ['On', 'On', antenas, freqs, direct_ZHS_inst]
#                   ]

# # type of simulation

# ZHAireS = True; RASPASS = False; upgoing = False 
# ZHAireS_control = [] if not ZHAireS else ZHAireS_control

# ########################### remote handling ############################

# jobID          = 'Simulator'
# remote_main    = '/home2/'
# user           = 'sergio.cabana/'
# remote_dir     = 'Pruebas_Running/Pruebas_submit/' # where to store files and future outputs
# exe            = ['aires/bin/ZHAireSRASPASS', 'SpecialPrimaries/RASPASSprimary', 'SpecialPrimaries/uprimary']
# program        = 'ZHAireSRASPASS'
# autorename     = True
# local_savepath = ''
# eco            = False


# server         = 'mastercr1.igfae.usc.es'
# node           = 'nodo014'
# username       = 'sergio.cabana'
# #########################################################################

# maz.simulator(task_names, basics, trajects, sim_controls, exports, extras, jobID, \
#               RASPASS = RASPASS, upgoing = upgoing, ZHAireS = ZHAireS, \
#               ZHAireS_control = ZHAireS_control, remote_main = remote_main, \
#               user = user, remote_dir = remote_dir, exe = exe, \
#               program = program, autorename = autorename, local_savepath = local_savepath, \
#               eco = eco, server = server, node = node, username = username)
