import os
import sys
sys.path.append('..')
from lib.simulation_objects import *
from lib.mppaper import *
import lib.mpsetup as mpsetup
import shutil
from scipy.integrate import solve_ivp
import matplotlib
import matplotlib.gridspec as gridspec
import json
import glob
import matplotlib.cm as cm
import matplotlib.colors as colors
import cmasher as cmr
import matplotlib.transforms as transforms
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
import matplotlib.font_manager as fm

from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import FancyArrow
from matplotlib.collections import PatchCollection

    
fontprops = fm.FontProperties(size=17)

dir_fig='../../sims_for_fig5' # directory with input files

dir_io='{inp}/t4_F8_vsF10'.format(inp=dir_fig)  # directory with input files

dir_in_tot='{inp}/data'.format(inp=dir_io) # directory with output plots
dir_in_frames='{inp}/frames'.format(inp=dir_in_tot) # directory with output plots
dir_out_plots_tot='{inp}/plots'.format(inp=dir_io) # directory with output plots
dir_out_plots_frames_int='{inp}/frames_intensity'.format(inp=dir_out_plots_tot) # directory with output plots


if not os.path.exists(dir_out_plots_tot):
    os.makedirs(dir_out_plots_tot) 
    

if not os.path.exists(dir_out_plots_frames_int):
    os.makedirs(dir_out_plots_frames_int) 
    
pars = {}
 
file_params='{inp}/parameters.json'.format(inp=dir_io)

with open(file_params,'r') as f_handle:
    pars = json.loads(f_handle.read())
    
print(pars)
print(type(pars))
print(pars['L'])



Ncols  = int(pars['L']/pars['Delta_x']) # number of grid rows and cols
Nrows  = int(pars['Ly']/pars['Delta_x']) # number of grid rows and cols
Nsites = int(Ncols*Nrows)

dx = dy = pars['Delta_x'] #lattice spacing, microns

discr_thr= (10.**8.)/(dx * dx *0.5) # discreteness treshold


L =pars['L'] 
model =pars['model'] 
print(model)


df = pd.read_excel('{inp}/../../growth_data/Scanner_OD.xlsx'.format(inp=dir_io), skiprows=1, usecols='A,B',  header = None,  engine='openpyxl', sheet_name="CFU Calibration")

pxtoCFUs=df.to_numpy(copy=True)[:8,:]

pxtoCFUs=pxtoCFUs[np.argsort(pxtoCFUs[:,1]),:]

print(pxtoCFUs)
print(pxtoCFUs[0,:])

# ~ f = interpolate.interp1d(pxtoCFUs[:,1], pxtoCFUs[:,0])
f = make_interp_spline(pxtoCFUs[:,1], pxtoCFUs[:,0], k=3)

############ PLOT FRAMES

times = []

bactprec=None    
intprec=None    
wavepos=None    
waveposprec=None    
frontpos=None    
frontposprec=None    
wavesum=0    


savedat=True
for file_in_frames in glob.glob("{inp}/experiment_state_time_*".format(inp=dir_in_frames)):
    filename, ext = os.path.splitext(file_in_frames)
    savedat=(ext=='.dat')
    print(ext, savedat)
    time=os.path.basename(filename).split('_')[-1]
    print(time)
    # ~ time=int(float(time))
    time=float(time)

    
    times.append(time)
    

times.sort()  

times =np.asarray(times)	
t_ref=[32]


for it, t in enumerate(times):
    
    
    if savedat:
        file_in_frames = "{inp}/experiment_state_time_{t}.dat".format(inp=dir_in_frames, t=t)
        data_space = np.loadtxt(file_in_frames)
    else:
        file_in_frames = "{inp}/experiment_state_time_{t}.npz".format(inp=dir_in_frames, t=t)
        with np.load(file_in_frames) as data:
             data_space = data['arr_0']
 
    print(t)



    x       = data_space[:,0]
    print (x)
    print (x.shape)
    y       = data_space[:,1]

    S   = data_space[:,3]
    E   = data_space[:,4]
    V = data_space[:,5]
    R=data_space[:,2]
    A=data_space[:,6]
    A2=data_space[:,7]
    
    B_tot = S+E
        
    S2 = np.zeros(x.shape)  # 
    Ei2 = np.zeros(E.shape)  # 

    if model == 'twobacts' and t>0:
        S2=data_space[:,8]
        Ei2=data_space[:,9]
        
        B_tot = B_tot + S2 + Ei2   
        
    xnew = B_tot
    intensities = f(xnew)
    intensities = np.array(intensities).astype(float)  
    
    print(intensities)
    print(intensities.dtype)
    
    
    # ~ xx = np.reshape(x, (Ncols, Ncols))/1000.
    xx = np.reshape(x, (Nrows, Ncols))/1000.
    print (xx.shape)
    # ~ sys.exit
    yy = np.reshape(y, xx.shape)/1000.
    R = np.reshape(R, xx.shape)
    S = np.reshape(S, xx.shape)
    E = np.reshape(E, xx.shape)
    V = np.reshape(V, xx.shape)
    A = np.reshape(A, xx.shape)
    A2 = np.reshape(A2, xx.shape)
    B_tot = np.reshape(B_tot, xx.shape)

    intensities = np.reshape(intensities, xx.shape)
 

    derivs_plots = derivatives_euler(pars, 'derivatives')
    
    newstate = derivs_plots(R,S,V,E, A,A2)  #solving the equation
    
    dRdt = newstate[0]
    dSdt = newstate[1]
    dVdt = newstate[2]
    dAdt = newstate[3]
    dA2dt = newstate[4]
    chemo_str = - newstate[5]
    diff_str = newstate[6]
    gr_str = newstate[7]
    inf_str = newstate[8]
    lys_str = newstate[9]
    attr_sense_field = newstate[10]
    
    grad_sensing=np.gradient(attr_sense_field, pars['Delta_x'])
    chemo_flux=B_tot*grad_sensing

         
    if it==0:
        # ~ np.unravel_index(np.argmax(a, axis=None), a.shape)
        id_first_bacts = np.unravel_index(np.argmax(B_tot, axis=None), B_tot.shape)[0]
        id_first_bacts_x = np.unravel_index(np.argmax(B_tot, axis=None), B_tot.shape)[1]
        id_first_virs  = np.unravel_index(np.argmax(V, axis=None), V.shape)[0] 
        # ~ id_first_virs  = np.argmax(V    , axis=0)
        
        print(id_first_bacts)
        print(np.unravel_index(np.argmax(B_tot, axis=None), B_tot.shape))
        # ~ sys.exit()


    
    
    if bactprec is None:
        wavesum = np.zeros_like(xx)
    else:
        
        # the front position, plotted in the figure, is calculated taking the difference between intensities at two consecutive saved times and isolationg the coordinates where the intensity decreases (therefore bacteria increase) and resources are still close to the initial concentration
        
        wave=intprec - intensities
        # ~ wavemask=(wave > 0) &(R>2.3e3)
        wavemask=(wave > 0) &(R>3.3e3)
        
        
        wave= np.where(wavemask,wave, 0.)
        
        wavesum = wavesum + wave
        
    bactprec=B_tot
    intprec=intensities
    waveposprec=wavepos
    frontposprec=frontpos

        
        
                
    if int(t) in t_ref:
        
        
        fig = plt.figure(figsize=(4,4))
        grid = gridspec.GridSpec(1, 1, left=0.03, right=0.97, top=0.95, bottom=0.03,
             wspace=0.4, hspace=0.35)
        labeled_axes = []
        ax = plt.Subplot(fig, grid[0, 0])
        fig.add_subplot(ax)
        labeled_axes.append(ax)
           
             
            
        heatmap = ax.pcolormesh(-xx, yy, intensities, cmap='gray',vmin=0, vmax=255, shading='auto', linewidth=0,rasterized=True) # faster than pcolor
        
        # PLOT initial conditions
        ax.plot(-(21000  -L/2.)/1000, (21000  -L/2.)/1000, marker='x', markersize=8, color='k', linestyle='')
        ax.plot(-(6000  -L/2.)/1000, (6000  -L/2.)/1000, marker='o',  markerfacecolor='none', markersize=8,   color='k', linestyle='')
        ax.plot((6000  -L/2.)/1000, (6000  -L/2.)/1000, marker='o',  markerfacecolor='none', markersize=8,   color='k', linestyle='')
        
        # PLOT fronts
        heatmap = ax.pcolormesh(-xx, yy,   wavesum , cmap='gray_r', alpha=0.5) # faster than pcolor
        
        ax.set_aspect('equal')

        scalebar = AnchoredSizeBar(ax.transData,
                               25, '25 mm', 'upper left', 
                               pad=0.2,
                               color='black',
                               frameon=False,
                               size_vertical=0.3, fontproperties=fontprops)
        ax.add_artist(scalebar)
      
        trans = transforms.blended_transform_factory(ax.transData, ax.get_xticklabels()[0].get_transform())

    #
        ax.tick_params(left = False, right = False , labelleft = False ,
                    labelbottom = False, bottom = False)
                    
        
        out_file='{out}/wave_timelapse_plus_lastframe.png'.format(out=dir_out_plots_tot, time=t)
        #    print out_file
        fig.savefig(out_file, dpi=300)
        fig.clf()
        plt.close('all')

    
    

