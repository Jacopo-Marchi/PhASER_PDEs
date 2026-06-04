

"""
Script Name: run_model.py
Description: runs the PhASERs model, takes the input/output directory as inline parameter
Author:      Jacopo Marchi
Date:        11/18/2025
Usage:       python run_model.py input_output_directory
Notes:       Please make sure that your shell calls the correct python distribution. The code expects to find a file in the indicated directory, specifying 5 key parameters. The optimal usage involves using a bash wrapper (see model_batch_runs.sh in this repository) to create the necessary directory and cycle through the main parameters, but this can be done by hand and this script can be called directly through python. The comments indicate portions of the code that need to be uncommented to produce specific figures. Importantly, select the correct parameters dictionary which specifies the parameters not set in the bash wrapper. This script requires the core functions defined in the simulations_object.py script in the library folder.

"""

import os
import sys
from lib.simulation_objects import *
import shutil
# ~ from scipy.integrate import solve_ivp
import json
import matplotlib.pyplot as plt
import datetime
import matplotlib.colors as colors


dir_io=sys.argv[1] # directory with input files

dir_out_tot='{inp}/data'.format(inp=dir_io) # directory with output data
dir_out_frames='{inp}/frames'.format(inp=dir_out_tot) # directory with output data
dir_out_plots_tot='{inp}/plots'.format(inp=dir_io) # directory with output plots
dir_out_plots_frames='{inp}/frames'.format(inp=dir_out_plots_tot) # directory with output plots

if not os.path.exists(dir_out_plots_tot):
    os.makedirs(dir_out_plots_tot) 
    
      

if not os.path.exists(dir_out_plots_frames):
    os.makedirs(dir_out_plots_frames) 
    
      

if not os.path.exists(dir_out_tot):
    os.makedirs(dir_out_tot) 
    
      

if not os.path.exists(dir_out_frames):
    os.makedirs(dir_out_frames) 
    
# LOAD PARAMETERS AND INITIALIZE SIMULATION
    
param_file='{inp}/params.txt'.format(inp=dir_io)


params = np.genfromtxt(param_file, dtype="|S10, |S10, |S10, |S10, |S10", names=['mf', 'Pc', 'gr_fac', 'phi', 'eta'])


print(params)
print(params.shape)




#PARAMETERS DICTIONARIES FOR THE MAIN SIMULATIONS PRESENTED IN THE PAPER, CORRESPONDING TO PARAMETERS LISTED IN THE SUPPLEMENTARY TEXT. UNCOMMENT THE DESIRED DICTIONARY. THE CORRESPONDING FIGURE/CONDITION IS SPECIFIED IN THE COMMENTS




#F9 and F10 25% agar, Fig 2

# ~ param_list={"model": "attractant", "eta": 5, "xi": 5500000.0, "eps": 1e-06, "phi": 2e-08, "D_B": 270000.0, "lambda": 70.0, "D_v": 10000.0, "r": 1.6, "k": 5000.0, "omega": 0.01, "am": 0.3, "ap": 100, "D_R": 3000000.0, "Delta_x": 100.0, "Delta_B": 1.0, "n": 5, "L": 200000.0, "Ly": 200000.0, "mot_factor": 1.0, "mu_1": 0.0, "mu_2": 3e-08, "k_attr": 1.5, "D_A": 3000000.0, "mu_1a": 0.0, "mu_2a": 1.4e-07, "k_attra": 3, "D_A2": 3000000.0, "xi2": 4000000.0, "Nstages": 10}# T7   F9 Fig 2A
# ~ param_list={"model": "attractant", "eta": 2, "xi": 4000000.0, "eps": 1e-06, "phi": 6e-08, "D_B": 180000.0, "lambda": 40.0, "D_v": 10000.0, "r": 1.6, "k": 5000.0, "omega": 0.01, "am": 0.3, "ap": 100, "D_R": 3000000.0, "Delta_x": 100.0, "Delta_B": 1.0, "n": 5, "L": 200000.0, "Ly": 200000.0, "mot_factor": 1.0, "mu_1": 0.0, "mu_2": 3e-08, "k_attr": 1.5, "D_A": 3000000.0, "mu_1a": 0.0, "mu_2a": 1.4e-07, "k_attra": 3, "D_A2": 3000000.0, "xi2": 4000000.0, "Nstages": 10}# T4 Vc1e7 , 6.7 mm/hr F10 Fig 4C, and F9 Fig2B



#F8 strain
# ~ param_list={"model": "attractant", "eta": 5, "xi": 6000000.0, "eps": 1e-06, "phi": 2e-08, "D_B": 180000.0, "lambda": 70.0, "D_v": 10000.0, "r": 1.3, "k": 7000.0, "omega": 0.01, "am": 0.3, "ap": 10, "D_R": 3000000.0, "Delta_x": 100.0, "Delta_B": 1.0, "n": 5, "L": 200000.0, "Ly": 200000.0, "mot_factor": 1.0, "mu_1": 0.0, "mu_2": 2e-08, "k_attr": 3.2, "D_A": 3000000.0, "mu_1a": 0.0, "mu_2a": 1.4e-07, "k_attra": 10, "D_A2": 3000000.0, "xi2": 6000000.0, "Nstages": 10} # F8 T7, Fig3
# ~ param_list={"model": "attractant", "eta": 5, "xi": 6000000.0, "eps": 1e-06, "phi": 2e-09, "D_B": 180000.0, "lambda": 50.0, "D_v": 10000.0, "r": 1.3, "k": 7000.0, "omega": 0.01, "am": 0.3, "ap": 10, "D_R": 3000000.0, "Delta_x": 100.0, "Delta_B": 1.0, "n": 5, "L": 200000.0, "Ly": 200000.0, "mot_factor": 1.0, "mu_1": 0.0, "mu_2": 2e-08, "k_attr": 3.2, "D_A": 3000000.0, "mu_1a": 0.0, "mu_2a": 1.4e-07, "k_attra": 10, "D_A2": 3000000.0, "xi2": 6000000.0, "Nstages": 10} # F8 T4, SI linear adsorp (Pc -1)

# ~ param_list={"model": "attractant", "eta": 2, "xi": 6000000.0, "eps": 1e-06, "phi": 6e-08, "D_B": 180000.0, "lambda": 40.0, "D_v": 10000.0, "r": 1.3, "k": 7000.0, "omega": 0.01, "am": 0.3, "ap": 10, "D_R": 3000000.0, "Delta_x": 100.0, "Delta_B": 1.0, "n": 5, "L": 200000.0, "Ly": 200000.0, "mot_factor": 1.0, "mu_1": 0.0, "mu_2": 2e-08, "k_attr": 3.2, "D_A": 3000000.0, "mu_1a": 0.0, "mu_2a": 1.4e-07, "k_attra": 10, "D_A2": 3000000.0, "xi2": 6000000.0, "Nstages": 10} # F8 T4  , Fig3

# ~ param_list={"model": "attractant", "eta": 2, "xi": 6000000.0, "eps": 1e-06, "phi": 6e-08, "D_B": 180000.0, "lambda": 40.0, "D_v": 10000.0, "r": 1.3, "k": 7000.0, "omega": 0.01, "am": 0.3, "ap": 10, "D_R": 3000000.0, "Delta_x": 100.0, "Delta_B": 1.0, "n": 5, "L": 120000.0, "Ly": 120000.0, "mot_factor": 1.0, "mu_1": 0.0, "mu_2": 2e-08, "k_attr": 3.2, "D_A": 3000000.0, "mu_1a": 0.0, "mu_2a": 1.4e-07, "k_attra": 10, "D_A2": 3000000.0, "xi2": 6000000.0, "Nstages": 10} # T4  smaller plate, use this for Fig 5 to merge colonies


# F8 AND F10 STRAINS DIFFERENT AGAR FIG 4
# ~ param_list={"model": "attractant", "eta": 5, "xi": 15200000.0, "eps": 1e-06, "phi": 2e-08, "D_B": 320000.0, "lambda": 70.0, "D_v": 10000.0, "r": 1.3, "k": 7000.0, "omega": 0.01, "am": 0.3, "ap": 10, "D_R": 3000000.0, "Delta_x": 100.0, "Delta_B": 1.0, "n": 5, "L": 200000.0, "Ly": 200000.0, "mot_factor": 1.0, "mu_1": 0.0, "mu_2": 2e-08, "k_attr": 3.2, "D_A": 3000000.0, "mu_1a": 0.0, "mu_2a": 1.4e-07, "k_attra": 10, "D_A2": 3000000.0, "xi2": 6000000.0, "Nstages": 10} # F8 T7, 0.2 agar, 8.3 mm/hr
# ~ param_list={"model": "attractant", "eta": 5, "xi": 2900000.0, "eps": 1e-06, "phi": 2e-08, "D_B": 100000.0, "lambda": 70.0, "D_v": 10000.0, "r": 1.3, "k": 7000.0, "omega": 0.01, "am": 0.3, "ap": 10, "D_R": 3000000.0, "Delta_x": 100.0, "Delta_B": 1.0, "n": 5, "L": 120000.0, "Ly": 120000.0, "mot_factor": 1.0, "mu_1": 0.0, "mu_2": 2e-08, "k_attr": 3.2, "D_A": 3000000.0, "mu_1a": 0.0, "mu_2a": 1.4e-07, "k_attra": 10, "D_A2": 3000000.0, "xi2": 6000000.0, "Nstages": 10} # F8 T7, 0.3 agar  3.5 mm/hr

# ~ param_list={"model": "attractant", "eta": 2, "xi": 15200000.0, "eps": 1e-06, "phi": 6e-08,  "D_B": 320000.0, "lambda": 40.0, "D_v": 10000.0, "r": 1.3, "k": 7000.0, "omega": 0.01, "am": 0.3, "ap": 10, "D_R": 3000000.0, "Delta_x": 100.0, "Delta_B": 1.0, "n": 5, "L": 200000.0, "Ly": 200000.0, "mot_factor": 1.0, "mu_1": 0.0, "mu_2": 2e-08, "k_attr": 3.2, "D_A": 3000000.0, "mu_1a": 0.0, "mu_2a": 1.4e-07, "k_attra": 10, "D_A2": 3000000.0, "xi2": 6000000.0, "Nstages": 10} # T4 0.2 agar,  8.3 mm/hr
# ~ param_list={"model": "attractant", "eta": 2, "xi": 2900000.0, "eps": 1e-06, "phi": 6e-08, "D_B": 100000.0, "lambda": 40.0, "D_v": 10000.0, "r": 1.3, "k": 7000.0, "omega": 0.01, "am": 0.3, "ap": 10, "D_R": 3000000.0, "Delta_x": 100.0, "Delta_B": 1.0, "n": 5, "L": 200000.0, "Ly": 200000.0, "mot_factor": 1.0, "mu_1": 0.0, "mu_2": 2e-08, "k_attr": 3.2, "D_A": 3000000.0, "mu_1a": 0.0, "mu_2a": 1.4e-07, "k_attra": 10, "D_A2": 3000000.0, "xi2": 6000000.0, "Nstages": 10} # F8 T4 0.3 agar,  3.5 mm/hr

# ~ param_list={"model": "attractant", "eta": 5, "xi": 4000000.0, "eps": 1e-06, "phi": 2e-08, "D_B": 180000.0, "lambda": 70.0, "D_v": 10000.0, "r": 1.6, "k": 5000.0, "omega": 0.01, "am": 0.3, "ap": 100, "D_R": 3000000.0, "Delta_x": 100.0, "Delta_B": 1.0, "n": 5, "L": 200000.0, "Ly": 200000.0, "mot_factor": 1.0, "mu_1": 0.0, "mu_2": 3e-08, "k_attr": 1.5, "D_A": 3000000.0, "mu_1a": 0.0, "mu_2a": 1.4e-07, "k_attra": 3, "D_A2": 3000000.0, "xi2": 4000000.0, "Nstages": 10}# F10 T7, 0.25 agar  6.7mm/hr , Fig 4C


#Different N_stages, named L in the manuscript. Set N_stages to the desired number
# ~ param_list={"model": "attractant", "eta": 5, "xi": 6000000.0, "eps": 1e-06, "phi": 2e-08, "D_B": 180000.0, "lambda": 70.0, "D_v": 10000.0, "r": 1.3, "k": 7000.0, "omega": 0.01, "am": 0.3, "ap": 10, "D_R": 3000000.0, "Delta_x": 100.0, "Delta_B": 1.0, "n": 5, "L": 600000.0, "Ly": 600000.0, "mot_factor": 1.0, "mu_1": 0.0, "mu_2": 2e-08, "k_attr": 3.2, "D_A": 3000000.0, "mu_1a": 0.0, "mu_2a": 1.4e-07, "k_attra": 10, "D_A2": 3000000.0, "xi2": 6000000.0, "Nstages": 1} # F8 T7 

# ~ param_list={"model": "attractant", "eta": 2, "xi": 6000000.0, "eps": 1e-06, "phi": 6e-08, "D_B": 180000.0, "lambda": 40.0, "D_v": 10000.0, "r": 1.3, "k": 7000.0, "omega": 0.01, "am": 0.3, "ap": 10, "D_R": 3000000.0, "Delta_x": 100.0, "Delta_B": 1.0, "n": 5, "L": 200000.0, "Ly": 200000.0, "mot_factor": 1.0, "mu_1": 0.0, "mu_2": 2e-08, "k_attr": 3.2, "D_A": 3000000.0, "mu_1a": 0.0, "mu_2a": 1.4e-07, "k_attra": 10, "D_A2": 3000000.0, "xi2": 6000000.0, "Nstages": 1} # F8 T4  




#F8 nutrient is attractant IMPORTANT, UNCOMMENT THE RESCALING OF am AND ap BELOW. NOTE THE DIFFERENT MODEL STRING
# ~ param_list={"model": "track_nutr", "eta": 5, "xi": 50000000.0, "eps": 1e-06, "phi": 2e-08, "D_B": 180000.0, "lambda": 70.0, "D_v": 10000.0, "r": 1.3, "k": 7000.0, "omega": 0.01, "am": 0.3, "ap": 10, "D_R": 3000000.0, "Delta_x": 100.0, "Delta_B": 1.0, "n": 5, "L": 200000.0, "Ly": 200000.0, "mot_factor": 1.0, "mu_1": 0.0, "mu_2": 2e-08, "k_attr": 3.2, "D_A": 3000000.0, "mu_1a": 0.0, "mu_2a": 1.4e-07, "k_attra": 10, "D_A2": 3000000.0, "xi2": 6000000.0, "Nstages": 10} # F8 T7 
# ~ param_list={"model": "track_nutr", "eta": 10, "xi": 50000000.0, "eps": 1e-06, "phi": 2e-08, "D_B": 180000.0, "lambda": 70.0, "D_v": 10000.0, "r": 1.3, "k": 7000.0, "omega": 0.01, "am": 0.3, "ap": 10, "D_R": 3000000.0, "Delta_x": 100.0, "Delta_B": 1.0, "n": 5, "L": 200000.0, "Ly": 200000.0, "mot_factor": 1.0, "mu_1": 0.0, "mu_2": 2e-08, "k_attr": 3.2, "D_A": 3000000.0, "mu_1a": 0.0, "mu_2a": 1.4e-07, "k_attra": 10, "D_A2": 3000000.0, "xi2": 6000000.0, "Nstages": 10} # F8 T7 
# ~ param_list["am"]    =  param_list["am"] * param_list["k"]/param_list["k_attra"] # Rescale so that sensitifity range is on the relevant order of magnitude as the chemical dynamics
# ~ param_list["ap"]    =  param_list["ap"] * param_list["k"]/param_list["k_attra"] # Rescale so that sensitifity range is on the relevant order of magnitude as the chemical dynamics




# SET/OVERRIDE A FEW REMAINING PARAMETERS

conv_LB=1. # # nutrients dilution in the initial LB preparation, 1 default
# ~ conv_LB = 0.50 #

param_list['mot_factor']=float(params['mf']) # motility factor for infected cells
param_list['Pc'] =float(params['Pc']) # phage adsorption saturation. If negative, a linear profile is employed (always the case for phage T7)
param_list['phi'] =float(params['phi']) # phage adsorption saturation. If negative, a linear profile is employed (always the case for phage T7)
param_list['eta'] =float(params['eta']) # phage adsorption saturation. If negative, a linear profile is employed (always the case for phage T7)
r_def = param_list['r'] # growth rate
param_list['r'] =float(params['gr_fac'])*r_def  # growth multiplied by modulation factor. 1 for default
param_list['eps'] =param_list['eps']/float(params['gr_fac'])  # yield multiplied by modulation factor. 1 for default. This assume that medium quality affects biomass conversion while chemicals uptake stays the same
param_list["xi2"]    = xi2 = param_list['xi'] # Ensures bacteria have an equal sensitivity to both attractant. Change if this is not desired


# DEFINE A FEW USEFUL VARIABLES
# ~ param_list['model'] = "twobacts" # UNCOMMENT TO GENERATE FIG 5E
model =param_list['model'] 

K =     param_list['k']

L  = param_list['L'] # petri dish size, in microns
Ly = param_list['Ly'] # petri dish size, in microns

dx = dy = Delta_x= param_list['Delta_x'] #lattice spacing, microns



Ncols  = int(param_list['L']/param_list['Delta_x']) # number of grid rows and cols
Nrows  = int(param_list['Ly']/param_list['Delta_x']) # number of grid rows and cols
Nsites = int(Ncols*Nrows)

Nstages = param_list['Nstages'] # number of infected stages

print(param_list["phi"])
# ~ sys.exit()
 
file_out='{inp}/parameters.json'.format(inp=dir_io)

with open(file_out,'w') as f_handle:
    json.dump(param_list, f_handle)
    

param_list_bact2={"model": "attractant", "eta": 2, "xi": 4000000.0, "eps": 1e-06, "phi": 6e-08, "D_B": 180000.0, "lambda": 40.0, "D_v": 10000.0, "r": 1.6, "k": 5000.0, "omega": 0.01, "am": 0.3, "ap": 100, "D_R": 3000000.0, "Delta_x": 100.0, "Delta_B": 1.0, "n": 5, "L": 200000.0, "Ly": 200000.0, "mot_factor": 1.0, "mu_1": 0.0, "mu_2": 3e-08, "k_attr": 1.5, "D_A": 3000000.0, "mu_1a": 0.0, "mu_2a": 1.4e-07, "k_attra": 3, "D_A2": 3000000.0, "xi2": 3500000.0, "Nstages": 10, "psi": 5.999999999999999e-06, "csi": 60.0}# Parameters for the second bacteria colony in figure 5E. Phage parameters are ignored here. This dictionary is only considered if mergefromscratch below is set to True and model is set to "twobacts" 

if model =="twobacts":
    
    file_out='{inp}/parameters_bact2.json'.format(inp=dir_io)

    with open(file_out,'w') as f_handle:
        json.dump(param_list_bact2, f_handle)
    


# LOAD CALIBRATION DATA TO CONVERT CFU TO PIXEL INTENSITIES

df = pd.read_excel('{inp}/../../growth_data/Scanner_OD.xlsx'.format(inp=dir_io), skiprows=1, usecols='A,B',  header = None,  engine='openpyxl', sheet_name="CFU Calibration")


pxtoCFUs=df.to_numpy(copy=True)[:8,:]

pxtoCFUs=pxtoCFUs[np.argsort(pxtoCFUs[:,1]),:]

print(pxtoCFUs)
print(pxtoCFUs[0,:])

f = make_interp_spline(pxtoCFUs[:,1], pxtoCFUs[:,0], k=3)

no_bactcontrol= np.amax(pxtoCFUs[:,0])

maxbacts_convert=np.amin(pxtoCFUs[:,0]) #min intensity

discr_thr= (10.**8.)/(Delta_x * Delta_x *0.5) # discreteness treshold for phage


# SET GRID AND INITIAL CONDITIONS

lindim=np.linspace(-param_list['L']/2, param_list['L']/2, Ncols, endpoint = False)
lindimy=np.linspace(-param_list['Ly']/2, param_list['Ly']/2, Nrows, endpoint = False)

xx, yy = np.meshgrid(lindim, lindimy)



# ~ A  = np.zeros(xx.shape) + 20.*conv_LB # attractant, IMPORTANT, UNCOMMENT FOR STRAINS F8 AND F10
A = np.zeros(xx.shape)  # attractant, uncomment and set this to 0 for strain F9
A2 = np.zeros(xx.shape) + 20.*conv_LB # attractant

print (A.shape)


R = np.zeros(xx.shape) + 3500*conv_LB # resources
Ei = np.zeros((Nstages, xx.shape[0], xx.shape[1]))  # infected cell at various stages


in_radius_bact= 2000
in_bact=8*10**8. # initial bacteria

 

in_x_bacts= int((6000 -L/2.)/dx)*dx # initial position of bacteria, always microns
in_y_bacts= int((6000 -L/2.)/dx)*dx# initial position of bacteria, always microns



in_radius_virs= 200
in_virs=10**9. # initial phage, PFU/ml
# ~ in_virs=0. # NO PHAGE
# ~ in_virs=5*10**5. # VARY TO TEST DUFFERENT INOCULI IN FIG 4

in_x_virs= int((21000  -L/2.)/dx)*dx # initial position of virs, always microns 
in_y_virs=  int((21000  -L/2.)/dx)*dx # initial position of virs, always microns


if param_list['Ly']<param_list['L']:
    in_y_virs= in_y_bacts =0
    in_x_virs=int((L*2/10.  -L/2.)/dx)*dx
    
V = in_virs*np.exp(-((xx- in_x_virs)**2. + (yy - in_y_virs) **2.)/(2*in_radius_virs **2.)) # initial phage distributed as a gaussian
S = in_bact*np.exp(-((xx- in_x_bacts)**2. + (yy - in_y_bacts) **2.)/(in_radius_bact **2.)) # initial bacteria distributed as a gaussian



S2 = np.zeros(xx.shape)  # 
Ei2 = np.zeros(Ei.shape)  # 

if model == 'twobacts':

    S2=S.copy()
    Ei2=Ei.copy()
    

    
# to merge colonies in Fig 5, set to True
# ~ mergefromscratch=True
mergefromscratch=False

if mergefromscratch:
    
    if model == 'twobacts': # laser two strains
        
        
        in_x_bacts2= - in_x_bacts  # initial position of virs, always microns
        in_y_bacts2= in_y_bacts

        
        S2 = in_bact*np.exp(-((xx- in_x_bacts2)**2. + (yy - in_y_bacts2) **2.)/(in_radius_bact **2.))
        
    else:
                
        # ~ in_x_virs= int((21000  -L/2.)/dx)*dx
        in_x_virs= int((25000  -L/2.)/dx)*dx
        # ~ in_x_bacts= int((5000 -L/2.)/dx)*dx # initial position of virs  always microns
        in_x_bacts= in_x_virs-(15000*np.sqrt(2)) # initial position of bacteria  always microns

    
        in_y_bacts= 0# initial position of bacteria, always microns
        in_y_virs= 0# initial position of virs, always microns
        
        
    
        xx_rot_tmp, yy_rot_tmp = rotate((0,0),xx, yy, np.deg2rad(45))
    
        xx_rot, yy_rot = rotate((70000  -L/2.,0),xx_rot_tmp, np.abs(yy_rot_tmp), np.deg2rad(60)) # half the angle between the PhASERs, set the desired angle (divided by two) HERE
        V = in_virs*np.exp(-((xx_rot- in_x_virs)**2. + (yy_rot - in_y_virs) **2.)/(2*in_radius_virs **2.))
    
        S = in_bact*np.exp(-((xx_rot- in_x_bacts)**2. + (yy_rot - in_y_bacts) **2.)/(in_radius_bact **2.))
        

#  to send PhASER to boundary in Supplementary Figure, set to True
# ~ crushagainstwall=True
crushagainstwall=False 

if crushagainstwall:
     
    in_x_virs= int((L/2. - 40000)/dx)*dx
    in_x_bacts= in_x_virs-(15000*np.sqrt(2)) # initial position of virs + 1cm shift, always microns


    in_y_bacts= 0# initial position of virs, always microns
    in_y_virs= 0# initial position of virs, always microns

    xx_rot, yy_rot = rotate((L/2.,0),xx, yy, np.deg2rad(-(90-40))) # angle between wall and PhASERs, set the desired one HERE
    
    V = in_virs*np.exp(-((xx_rot- in_x_virs)**2. + (yy_rot - in_y_virs) **2.)/(2*in_radius_virs **2.))
    S = in_bact*np.exp(-((xx_rot- in_x_bacts)**2. + (yy_rot - in_y_bacts) **2.)/(in_radius_bact **2.))
    





elapsed = 0.
elapsed_print = -100.
elapsed_check = -100. # with negative values densities are checked at every iteration. Set higher to speed up code, but paying attention not to introduce numerical errors.
totalTime = 200

init_realtime=datetime.datetime.now()


dt = 5e-3  #starting dt of the simulation.

neg_R = 0 # total number of times resources went below allowed limit. Useful to refine time stepping
neg_P = 0
neg_S = 0
neg_E1 = 0
neg_E2 = 0
neg_A = 0
neg_Vs = 0

store_timeStep = [dt]
store_elapsed = [elapsed]

store_steps_all = [0]


store_neg_A =  [0]
store_neg_R =  [0] # stores how many time iterations had resources below allowed limit. Useful to refine time stepping
store_neg_P =  [0]
store_neg_S =  [0]
store_neg_E1 = [0]
store_neg_E2 = [0]

steps = 0
steps_all = 0


            
Ei_flat = np.reshape(np.ravel(Ei), (Nstages, Nsites))

Ei_all= np.sum(Ei_flat, axis=0) # I want to save all infected together

Ei2_flat = np.reshape(np.ravel(Ei2), (Nstages, Nsites))

Ei2_all= np.sum(Ei2_flat, axis=0) # I want to save all infected together

print(Ei.shape)
print(Ei_all.shape)

xs = np.ravel(xx) # already in microns
ys = np.ravel(yy) # 

print(S.shape)

# store output 

data_fin = np.array([xs, ys,np.ravel(R), np.ravel(S), Ei_all, np.ravel(V), np.ravel(A) , np.ravel(A2), np.ravel(S2), Ei2_all ])
savedat=False # SET TO TRUE TO SAVE OUTPUT IN READABLE TEXT FILES. Much heavier, but useful during debugging. Otherwise output is compressed in binary npz files 
if savedat:
    file_out='{inp}/experiment_state_time_{t}.dat'.format(inp=dir_out_frames, t=elapsed)
    
    with open(file_out,'a') as f_handle:
        np.savetxt(f_handle, data_fin.T, fmt='%40.15f',  header="x \t y \t R \t S \t Ei_all \t V")
        f_handle.write("\n")
else:
    file_out='{inp}/experiment_state_time_{t}.npz'.format(inp=dir_out_frames, t=elapsed)
    
    np.savez_compressed(file_out, data_fin.T)


B=S+np.sum(Ei, axis =0 ) +np.sum(Ei2, axis =0 ) + S2


xnew = B.ravel()
intensities = f(xnew)
intensities = np.array(intensities).astype(float)  
intensities = np.reshape(intensities, xx.shape)

epsilon_plot = 1000.

# PLOT INITIAL CONDITIONS

fig = plt.figure(figsize=(8,8), constrained_layout=False)
grid = fig.add_gridspec(nrows=3,ncols=2, left=0.1, right=0.95, top=0.97, bottom=0.04,
         wspace=0.4, hspace=0.2)
axint = fig.add_subplot(grid[0:-1, :])


    
heatmap = axint.pcolormesh(xx/epsilon_plot, yy/epsilon_plot, intensities, cmap='gray',vmin=0, vmax=255,shading='auto', linewidth=0,rasterized=True) # faster than pcolor
axint.set_xlabel(r'x (mm)')
axint.set_ylabel(r'y (mm)')

axint.axis('scaled')

axint.set_title("Time: " +str(int(elapsed))+" hours", fontsize=16)

    
ax0 = fig.add_subplot(grid[-1, 1])
ax1 = fig.add_subplot(grid[-1, 0])

v1 = ax0.pcolormesh(xx/epsilon_plot,yy/epsilon_plot, V,norm=colors.LogNorm(vmin=discr_thr, vmax=1e11),shading='auto',cmap='Oranges', linewidth=0,rasterized=True)  ##
v2 = ax1.pcolormesh(xx/epsilon_plot,yy/epsilon_plot, B, norm=colors.LogNorm(vmin=discr_thr, vmax=1e10),shading='auto',cmap = 'PuBu', linewidth=0,rasterized=True) 

ax0.set_xlabel(r'x (mm)')
ax0.set_ylabel(r'y (mm)')

ax1.set_xlabel(r'x (mm)')
ax1.set_ylabel(r'y (mm)')

ax0.axis('scaled')
ax1.axis('scaled')



cbar0 = fig.colorbar(heatmap,ax=axint)
cbar0.ax.get_yaxis().labelpad = 15
cbar0.ax.set_ylabel('pixel intensity', rotation=270)
    


cbar1 = fig.colorbar(v1,ax=ax0)
cbar2 = fig.colorbar(v2,ax=ax1)

cbar1.ax.get_yaxis().labelpad = 15
cbar1.ax.set_ylabel(r"Phage density $(\rm{PFU} / mL)$", rotation=270)
cbar2.ax.get_yaxis().labelpad = 15
cbar2.ax.set_ylabel(r"Bacteria density $(\rm{CFU} / mL)$", rotation=270)

fig.savefig(dir_out_plots_frames+"/"+"frame_t_{t}.pdf".format(t=elapsed))



fig.clf()
plt.close('all')

R_new =  R
S_new = S
V_new = V
Ei_new =Ei
A_new = A
A2_new =A2
S2_new =S2
Ei2_new =Ei2


# DEFINE FUNCTION TO ADVANCE SIMULATION TO THE NEXT STEP
advance = derivatives_euler(param_list, model)

if model =="twobacts":
    advance = derivatives_euler(param_list, model, param_list_bact2)
    

checkfreq=0.0 #hours

# INTEGRATE PDES UNTIL THE MAXIMUM TIME IS REACHED OR THE TIME STEP GOES BELOW 1E-6
while elapsed < totalTime:
    
    steps_all= steps_all+1
    
    # ~ print elapsed, dt
    
    if dt < 1e-6 :
        print("dt too small")
        
        
        now = datetime.datetime.now()
        print("Current real time : ")
        print(now.strftime("%H:%M:%S"))

        print("Real Time Elapsed ", (now - init_realtime).total_seconds() /60., " min")
        print("Current dt: ", dt)
        print("Simulation Time Elapsed: ", elapsed, " ", steps)
        print("Bacteria max value: ",np.max(S),"  Bacteria min value: ", np.min(S))
        print("Phage max value: ",np.max(V),"  Phage min value: ", np.min(V))
        print("Resources max value: ",np.max(R),"  Resources min value: ", np.min(R))
        print("E1 max and min values: ",np.max(Ei), np.min(Ei))
        print("Negative S count: ", neg_S)
        print("Negative P count: ", neg_P)
        print("Negative R count: ", neg_R)
        print("Negative Ei count: ", neg_E1)
        if model == 'attractant' or model == 'twobacts':
            print("Negative A count: ", neg_A)
        print("\n")
        break
    

    dt_until = totalTime - elapsed   #checking how much time is left for the simulation
    dt_save = dt
    
    if dt_until < dt:
        dt = dt_until
        
    

    if model == 'attractant':
        newstate = advance(dt, R,S,V,Ei, A,A2)  #advance the state variables to the next step
        
        R_new = newstate[0]
        S_new = newstate[1]
        V_new = newstate[2]
        Ei_new = newstate[3]
        A_new = newstate[4]
        A2_new = newstate[5]
    elif model == 'twobacts':
        newstate = advance(dt, R,S,V,Ei, A,A2, S2, Ei2)  #advance the state variables to the next step
        
        R_new = newstate[0]
        S_new = newstate[1]
        V_new = newstate[2]
        Ei_new = newstate[3]
        A_new  = newstate[4]
        A2_new = newstate[5]
        S2_new = newstate[6]
        Ei2_new = newstate[7]
        
    elif model == 'track_nutr':
        newstate = advance(dt, R,S,V,Ei)  #advance the state variables to the next step
        
        R_new = newstate[0]
        S_new = newstate[1]
        V_new = newstate[2]
        Ei_new = newstate[3]
    

# COUNT VARIABLES BELOW NEGATIVE THRESHOLDS, TO DETERMINE WHETHER TO LOWER TIME STEP     
    c1 = 0 
    c2 = 0 # obsolete
    c3 = 0 
    c4 = 0 # obsolete
    c5 = 0.  
    
    if (elapsed - elapsed_check)>=checkfreq:
        
        c1 = np.count_nonzero(R_new   < -K/1000)  #check if there are any negative numbers in R

        c3 = np.count_nonzero(V_new < -discr_thr/10) #check if there are any negative numbers in V
        if model == 'attractant' or model == 'twobacts':
            c5 = np.count_nonzero(A_new   < - 3*param_list["am"] ) + np.count_nonzero(A2_new   < -3*param_list["am"])  # #check if there are any negative numbers in A
        
        elapsed_check= elapsed
    
    if c1 != 0:
        neg_R += 1
    if c2 != 0:
        neg_S += 1
    if c3 != 0:
        neg_P += 1
    if c4 != 0:
        neg_E1 += 1
    
    if c5 != 0:
        neg_A += 1
    
    
    
    if c1 != 0 or c2 != 0 or c3 != 0 or c4 != 0 or c5 != 0: # IF ANY VARIABLE GOES BELOW ALLOWED LIMIT, DISREGARD LAST STEP AND REDUCE TIME STEP
        
        dt = dt / 1.5         #reduce dt to  improve the numerical integration, 
        
        checkfreq=0 #hours, check every step

    

    else: # IF NO VARIABLE GOES TOO FAR BELOW 0 (SEE ABOVE), RETAIN THE LAST STEP AND STORE THE OUTPUT
        
        checkfreq=0.0 #hours
        
        # set negative densities to 0. The check above guatrrantees that the time step is adapted so that densities are never extrapolated to significantly negative values
        R_new[R_new < 0 ] = 0.
        Ei_new[Ei_new < 0. ] = 0.
        S_new[S_new < 0. ] = 0.
        if model == 'attractant' or model == 'twobacts':
            A_new[A_new < 0. ] = 0.
            A2_new[A2_new < 0. ] = 0.
        
        if model == 'twobacts':
            S2_new[S2_new < 0. ] = 0.
            Ei2_new[Ei2_new < 0. ] = 0.
        
        V_new[V_new < 0. ] = 0.
        
        
    # copy  new states to current, and update time
        np.copyto(R , R_new  )
        np.copyto(S , S_new  )
        np.copyto(V , V_new  )
        np.copyto(Ei, Ei_new ) 
        if model == 'attractant' or model == 'twobacts':
            np.copyto(A,A_new)  
            np.copyto(A2 ,A2_new)  
        if model == 'twobacts':
            np.copyto(S2,S2_new)  
            np.copyto(Ei2,Ei2_new)  
        elapsed+=dt
        
        if (elapsed - elapsed_print)>=2: # save output every 2 hours
    
            steps = steps + 1
    
            store_elapsed.append(elapsed)
            store_timeStep.append(dt)
            
            
            store_steps_all.append(steps_all)
            store_neg_R.append(neg_R) # stores how many time iterations had resources below allowed limit. Useful to refine time stepping
            store_neg_P.append(neg_P )
            store_neg_S.append(neg_S )
            store_neg_E1.append(neg_E1)
            store_neg_A.append(neg_A)
            
            now = datetime.datetime.now()
            print("Current real time : ")
            print(now.strftime("%H:%M:%S"))
    
            # print elapsed time in microseconds
            print("Real Time Elapsed ", (now - init_realtime).total_seconds() /60., " min")
            print("Current dt: ", dt)
            print("Simulation Time Elapsed: ", elapsed, " ", steps)
            print("Bacteria max value: ",np.max(S),"  Bacteria min value: ", np.min(S))
            print("Phage max value: ",np.max(V),"  Phage min value: ", np.min(V))
            print("Resources max value: ",np.max(R),"  Resources min value: ", np.min(R))
            print("E1 max and min values: ",np.max(Ei), np.min(Ei))
            print("Negative S count: ", neg_S)
            print("Negative P count: ", neg_P)
            print("Negative R count: ", neg_R)
            print("Negative Ei count: ", neg_E1)
            if model == 'attractant' or model == 'twobacts':
                print("Negative A count: ", neg_A)
            
            print("\n")
            
            elapsed_print= elapsed
            
            Ei_flat = np.reshape(np.ravel(Ei), (Nstages, Nsites))
        
            Ei_all= np.sum(Ei_flat, axis=0) # I want to save all infected together
            
            Ei2_flat = np.reshape(np.ravel(Ei2), (Nstages, Nsites))
        
            Ei2_all= np.sum(Ei2_flat, axis=0) # I want to save all infected together
            
            print(Ei.shape)
            print(Ei_all.shape)
            
            xs = np.ravel(xx) # already in microns
            ys = np.ravel(yy) # 
            
            print(S.shape)
            
            # store output 
            data_fin = np.array([xs, ys,np.ravel(R), np.ravel(S), Ei_all, np.ravel(V), np.ravel(A) , np.ravel(A2), np.ravel(S2), Ei2_all ])
            if savedat:
                file_out='{inp}/experiment_state_time_{t}.dat'.format(inp=dir_out_frames, t=elapsed)
                
                with open(file_out,'a') as f_handle:
                    np.savetxt(f_handle, data_fin.T, fmt='%40.15f',  header="x \t y \t R \t S \t Ei_all \t V")
                    f_handle.write("\n")
            else:
                file_out='{inp}/experiment_state_time_{t}.npz'.format(inp=dir_out_frames, t=elapsed)
                
                np.savez_compressed(file_out, data_fin.T)

            
            #PLOT
            print("PLOTTING")
            
            B=S+np.sum(Ei, axis =0 )
            if model == 'twobacts':
                B=B+S2+np.sum(Ei2, axis =0 )

            
            xnew = B.ravel()
            intensities = f(xnew)
            intensities = np.array(intensities).astype(float)  
            intensities = np.reshape(intensities, xx.shape)
            
            epsilon_plot = 1000.
            
            
            fig = plt.figure(figsize=(8,8), constrained_layout=False)
            grid = fig.add_gridspec(nrows=3,ncols=2, left=0.1, right=0.95, top=0.97, bottom=0.04,
                     wspace=0.4, hspace=0.2)
            axint = fig.add_subplot(grid[0:-1, :])
            
            
                
            heatmap = axint.pcolormesh(xx/epsilon_plot, yy/epsilon_plot, intensities, cmap='gray',vmin=0, vmax=255,shading='auto', linewidth=0,rasterized=True) # faster than pcolor
            axint.set_xlabel(r'x (mm)')
            axint.set_ylabel(r'y (mm)')
            
            cbar1 = fig.colorbar(heatmap,ax=axint)
            
            cbar1.ax.get_yaxis().labelpad = 15
            cbar1.ax.set_ylabel('pixel intensity', rotation=270)
                
            
            axint.set_aspect('equal', adjustable='box')
            
            axint.set_title("Time: " +str(int(elapsed))+" hours", fontsize=16)
            ax0 = fig.add_subplot(grid[-1, 1])
            ax1 = fig.add_subplot(grid[-1, 0])
            
            v1 = ax0.pcolormesh(xx/epsilon_plot,yy/epsilon_plot, np.where(V>discr_thr, V,0) ,norm=colors.LogNorm(vmin=discr_thr, vmax=1e11),shading='auto',cmap='Oranges', linewidth=0,rasterized=True)  ##plotting with log scale
            v2 = ax1.pcolormesh(xx/epsilon_plot,yy/epsilon_plot, np.where(B>discr_thr, B,0), norm=colors.LogNorm(vmin=discr_thr, vmax=1e10),shading='auto',cmap = 'PuBu', linewidth=0,rasterized=True) 
            
            ax0.set_xlabel(r'x (mm)')
            ax0.set_ylabel(r'y (mm)')
            
            ax1.set_xlabel(r'x (mm)')
            ax1.set_ylabel(r'y (mm)')
            
            
            
            cbar1 = fig.colorbar(v1,ax=ax0)
            cbar2 = fig.colorbar(v2,ax=ax1)
            cbar1.ax.get_yaxis().labelpad = 15
            cbar1.ax.set_ylabel(r"Phage density $(\rm{PFU} / mL)$", rotation=270)
            cbar2.ax.get_yaxis().labelpad = 15
            cbar2.ax.set_ylabel(r"Bacteria density $(\rm{CFU} / mL)$", rotation=270)
            ax0.axis('scaled')
            ax1.axis('scaled')
            fig.savefig(dir_out_plots_frames+"/"+"frame_t_{t}.pdf".format(t=elapsed))
            
            fig.clf()
            plt.close('all')
                
        dt = dt_save
        if dt < 5e-4: 
            dt = dt * 1.3    #increase dt by a factor 1.3, up to a maximum of 5e-4


    

""" 
 Creating a new directory for the current simulation and it's specific parameters.
 The directory name will have some of the parameters (do we need to add all?)
 
 Saving the different 'store_X' arrays into different text files.
 
 Also saving a json file with the list of parameters for a better reference.
 
"""



prefix = "store"
names = {0: "_timeStep", 1: "_elapsed", 2: "_totsteps", 3: "_neg_R", 4: "_neg_P", 5: "_neg_S", 6: "_neg_Ei", 7: "_neg_A"}
arr = {0: store_timeStep, 1: store_elapsed, 2: store_steps_all, 3: store_neg_R, 4: store_neg_P, 5: store_neg_S, 6: store_neg_E1, 7: store_neg_A}

for i in range(8):
    np.savetxt(dir_out_tot +"/" + prefix + names[i] +".txt", arr[i])
                
                

