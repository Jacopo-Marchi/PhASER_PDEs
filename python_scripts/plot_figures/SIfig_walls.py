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


def plot_colourline(x,y,c,ax):
    col = cm.jet((c-np.min(c))/(np.max(c)-np.min(c)))
    # ~ ax = plt.gca()
    for i in np.arange(len(x)-1):
        ax.plot([x[i],x[i+1]], [y[i],y[i+1]], c=col[i])
    im = ax.scatter(x, y, c=c, s=0, cmap=cm.jet)
    return im

def plot_timelapse_intensity(xs,ys,zs,ax,idx_lapse,rightside=None,leftside=None,yoffset=None,rescdim=1, mfs_tit=None):
    v = ax.pcolormesh(xs, ys, zs, cmap='gray',vmin=0, vmax=255,shading='auto', linewidth=0,rasterized=True)
    # ~ axs_timelapse[i_lapse].set_title(str(take_elapsed[index[i]])+" hours", fontsize=17)
    ax.set_aspect('equal')
    # ~ axs_timelapse[i_lapse].tick_params(axis='x', labelsize=14)
    # ~ axs_timelapse[i_lapse].tick_params(axis='y', labelsize=14)
    
    

    if idx_lapse==0:
        scalebar = AnchoredSizeBar(ax.transData,
                               25, '', 'lower left', 
                               pad=0.2,
                               color='black',
                               frameon=False,
                               size_vertical=0.3, fontproperties=fontprops)
        ax.add_artist(scalebar)
    # ~ ax.legend(frameon=False, columnspacing=0.5, handletextpad=0.2,
      # ~ loc='upper right', bbox_to_anchor=(2.5, 1.2))
      
    trans = transforms.blended_transform_factory(ax.transData, ax.get_xticklabels()[0].get_transform())

    #
    ax.tick_params(left = False, right = False , labelleft = False ,
                labelbottom = False, bottom = False)
                
        
    if rightside is not None:
        width_inset=rightside-leftside
                    
        rect = patches.Rectangle((leftside, leftside+yoffset), width_inset, width_inset, linewidth=2, edgecolor='k', facecolor='none')

        # Add the patch to the Axes
        ax.add_patch(rect)
        
    return v

    
fontprops = fm.FontProperties(size=12)

dir_fig='../../sims_for_fig4' # directory with input files
dir_io='{inp}/crushwall_40deg'.format(inp=dir_fig)  # directory with input files

dir_in_tot='{inp}/data'.format(inp=dir_io) # directory with output plots
dir_in_frames='{inp}/frames'.format(inp=dir_in_tot) # directory with output plots
dir_out_plots_tot='{inp}/plots'.format(inp=dir_io) # directory with output plots


if not os.path.exists(dir_out_plots_tot):
    os.makedirs(dir_out_plots_tot) 
    



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

no_bactcontrol= np.amax(pxtoCFUs[:,0])

maxbacts_convert=np.amin(pxtoCFUs[:,0]) #min intensity
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
    
        


xx               =None
yy               =None
R                =None
S                =None
E                =None
V                =None
A                =None
A2               =None
BR               =None
B_tot            =None
intensities      =None
dRdt             =None
dSdt             =None
dVdt             =None
dAdt             =None
dA2dt            =None
chemo_str        =None
diff_str         =None
gr_str           =None
inf_str          =None
lys_str          =None
attr_sense_field =None
grad_sensing     =None
chemo_flux       =None


times.sort()  

times =np.asarray(times)	
t_ref=[-1]
t_int=[36]

    




for it, t in enumerate(times):
    
    if int(t) in t_int + t_ref :

        
        if savedat:
            file_in_frames = "{inp}/experiment_state_time_{t}.dat".format(inp=dir_in_frames, t=t)
            # ~ filename, _ = os.path.splitext(file_in_frames)
            # ~ time=os.path.basename(filename).split('_')[-1]
            # ~ time=int(float(time))
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
        BR = np.zeros(x.shape)  # 
    
        if model == 'resistance' and t>0:
            BR=data_space[:,8]
            
            B_tot = B_tot + BR    
    
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
        BR = np.reshape(BR, xx.shape)
        B_tot = np.reshape(B_tot, xx.shape)
    
        intensities = np.reshape(intensities, xx.shape)
    
        B_tot = S+E
        BR = np.zeros(xx.shape)  # 
    
        if model == 'resistance' and t>0:
            BR=data_space[:,8]
            BR = np.reshape(BR, xx.shape)
            
            B_tot = B_tot + BR
            
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
            
            
            
            # ~ wavemask=(intensities- intprec > 0) &(V<1e7)
            # ~ wave=intensities - intprec
            
            wave=intprec - intensities
            # ~ wavemask=(wave > 0) &(R>2.3e3)
            wavemask=(wave > 0) &(R>3.3e3)
            
            
            wave= np.where(wavemask,wave, 0.)
            
            wavesum = wavesum + wave
            
            
            
        bactprec=B_tot
        intprec=intensities
        waveposprec=wavepos
        frontposprec=frontpos
    

    if int(t) in t_int:
        
        fig = plt.figure(figsize=(4,4))
        grid = gridspec.GridSpec(1, 1, left=0.03, right=0.97, top=0.95, bottom=0.03,
             wspace=0.4, hspace=0.35)
        labeled_axes = []
        ax = plt.Subplot(fig, grid[0, 0])
        fig.add_subplot(ax)
        

        v= plot_timelapse_intensity(xx,yy,intensities,ax,0) #,rescdim=10/6
        
        in_x_virs= int((L/2. - 40000)/dx)*dx
        in_x_bacts= in_x_virs-(15000*np.sqrt(2)) # initial position of virs + 1cm shift, always microns
    
    
        in_y_bacts= 0# initial position of virs, always microns
        in_y_virs= 0# initial position of virs, always microns
    
        xx_rot, yy_rot = rotate((L/2.,0),xx, yy, np.deg2rad(-(90-40)))

        
        xx_rot_vir, yy_rot_vir = rotate((L/2.,0),in_x_virs, in_y_virs, np.deg2rad((90-40)))
        xx_rot_bacts, yy_rot_bacts = rotate((L/2.,0),in_x_bacts, in_y_bacts,  np.deg2rad((90-40)))
        
              
        ax.plot(xx_rot_vir/1000, yy_rot_vir/1000, marker='x', markersize=6, color='k', linestyle='')
        ax.plot(xx_rot_bacts/1000, yy_rot_bacts/1000, marker='o',  markerfacecolor='none', markersize=6,   color='k', linestyle='')
        
  
        
        out_file='{out}/last_intensity.png'.format(out=dir_out_plots_tot)
        #    print out_file
        fig.savefig(out_file, dpi=300)
        fig.clf()
        plt.close('all')
   
    
 
    
dir_io='{inp}/crushwall_60deg'.format(inp=dir_fig)  # directory with input files

dir_in_tot='{inp}/data'.format(inp=dir_io) # directory with output plots
dir_in_frames='{inp}/frames'.format(inp=dir_in_tot) # directory with output plots
dir_out_plots_tot='{inp}/plots'.format(inp=dir_io) # directory with output plots


if not os.path.exists(dir_out_plots_tot):
    os.makedirs(dir_out_plots_tot) 
    



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

no_bactcontrol= np.amax(pxtoCFUs[:,0])

maxbacts_convert=np.amin(pxtoCFUs[:,0]) #min intensity
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
t_ref=[-1]
t_int=[34]
    
    


for it, t in enumerate(times):
    
    if int(t) in t_int + t_ref:
        
        if savedat:
            file_in_frames = "{inp}/experiment_state_time_{t}.dat".format(inp=dir_in_frames, t=t)
            # ~ filename, _ = os.path.splitext(file_in_frames)
            # ~ time=os.path.basename(filename).split('_')[-1]
            # ~ time=int(float(time))
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
        BR = np.zeros(x.shape)  # 
    
        if model == 'resistance' and t>0:
            BR=data_space[:,8]
            
            B_tot = B_tot + BR    
    
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
        BR = np.reshape(BR, xx.shape)
        B_tot = np.reshape(B_tot, xx.shape)
    
        intensities = np.reshape(intensities, xx.shape)
    
        B_tot = S+E
        BR = np.zeros(xx.shape)  # 
    
        if model == 'resistance' and t>0:
            BR=data_space[:,8]
            BR = np.reshape(BR, xx.shape)
            
            B_tot = B_tot + BR
            
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
            
            
            
            # ~ wavemask=(intensities- intprec > 0) &(V<1e7)
            # ~ wave=intensities - intprec
            
            wave=intprec - intensities
            # ~ wavemask=(wave > 0) &(R>2.3e3)
            wavemask=(wave > 0) &(R>3.3e3)
            
            
            wave= np.where(wavemask,wave, 0.)
            
            wavesum = wavesum + wave
            
            
            
        bactprec=B_tot
        intprec=intensities
        waveposprec=wavepos
        frontposprec=frontpos


    if int(t) in t_int:
        
        fig = plt.figure(figsize=(4,4))
        grid = gridspec.GridSpec(1, 1, left=0.03, right=0.97, top=0.95, bottom=0.03,
             wspace=0.4, hspace=0.35)
        labeled_axes = []
        ax = plt.Subplot(fig, grid[0, 0])
        fig.add_subplot(ax)
        

        v= plot_timelapse_intensity(xx,yy,intensities,ax,0) #,rescdim=10/6
        
        in_x_virs= int((L/2. - 40000)/dx)*dx
        in_x_bacts= in_x_virs-(15000*np.sqrt(2)) # initial position of virs + 1cm shift, always microns
    
    
        in_y_bacts= 0# initial position of virs, always microns
        in_y_virs= 0# initial position of virs, always microns
    

        
        # ~ in_y_bacts= 10000# initial position of virs, always microns
        # ~ in_y_virs= 10000# initial position of virs, always microns
    

        # ~ xx_rot, yy_rot = rotate((76000  -L/2.,0),xx_rot_tmp, np.abs(yy_rot_tmp), np.deg2rad(22.5))
        xx_rot_vir, yy_rot_vir = rotate((L/2.,0),in_x_virs, in_y_virs, np.deg2rad((90-60)))
        xx_rot_bacts, yy_rot_bacts = rotate((L/2.,0),in_x_bacts, in_y_bacts,  np.deg2rad((90-60)))
        
              
        ax.plot(xx_rot_vir/1000, yy_rot_vir/1000, marker='x', markersize=6, color='k', linestyle='')
        ax.plot(xx_rot_bacts/1000, yy_rot_bacts/1000, marker='o',  markerfacecolor='none', markersize=6,   color='k', linestyle='')
        
  
        
        out_file='{out}/last_intensity.png'.format(out=dir_out_plots_tot)
        #    print out_file
        fig.savefig(out_file, dpi=300)
        fig.clf()
        plt.close('all')
        
    

