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



def plot_timelapse_intensity(xs,ys,zs,ax,idx_lapse,rightside=None,leftside=None,yoffset=None,rescdim=1, mfs_tit=None):
    v = ax.pcolormesh(xs, ys, zs, cmap='gray',vmin=0, vmax=255,shading='auto', linewidth=0,rasterized=True)
    # ~ axs_timelapse[i_lapse].set_title(str(take_elapsed[index[i]])+" hours", fontsize=17)
    ax.set_aspect('equal')
    

    if idx_lapse==0:
        scalebar = AnchoredSizeBar(ax.transData,
                               25, '25 mm', 'lower left', 
                               pad=0.2,
                               color='black',
                               frameon=False,
                               size_vertical=0.3, fontproperties=fontprops)
        ax.add_artist(scalebar)
    # ~ ax.legend(frameon=False, columnspacing=0.5, handletextpad=0.2,
      # ~ loc='upper right', bbox_to_anchor=(2.5, 1.2))
      
    trans = transforms.blended_transform_factory(ax.transData, ax.get_xticklabels()[0].get_transform())
    if mfs_tit is not None:
        if '\n' not in mfs_tit:
            ax.annotate(mfs_tit, xy=(-90/rescdim, 85/rescdim), xycoords='data',xytext=(0, 5), textcoords='offset pixels', color='k', fontsize=16)
        else:
            ax.annotate(mfs_tit, xy=(-90/rescdim, 75/rescdim), xycoords='data',xytext=(0, 5), textcoords='offset pixels', color='k', fontsize=16)
    #
    ax.tick_params(left = False, right = False , labelleft = False ,
                labelbottom = False, bottom = False)
                
        
    if rightside is not None:
        width_inset=rightside-leftside
                    
        rect = patches.Rectangle((leftside, leftside+yoffset), width_inset, width_inset, linewidth=2, edgecolor='k', facecolor='none')

        # Add the patch to the Axes
        ax.add_patch(rect)
        
    return v

def plot_quiver(xs,ys,zs,ax,idx_lapse):

    
    step=6
    heatmap_int = ax.pcolormesh(xs, ys, zs, cmap='gray',vmin=0, vmax=255,shading='auto', linewidth=0,rasterized=True) # faster than pcolor
    
    cutvir=np.amax(V)/20.
    V_cut= np.where(V>cutvir,V,np.nan)
    

    heatmap_V = ax.pcolormesh(xs, ys, V_cut, cmap='Reds',vmin=1e6, vmax=1e10, alpha=0.5) # faster than pcolor, norm=colors.LogNorm()
  
    # ~ cutbact=min(discr_thr, np.amax(E)/20.)
    cutbact=np.amax(E)/20.
    E_cut= np.where(E>cutbact,E,np.nan)

    heatmap_E = ax.pcolormesh(xx, yy, E_cut, cmap='PuBu',vmin=1e6, vmax=1e9, alpha=0.5) # faster than pcolor, norm=colors.LogNorm()
        
    values_cont=[pars['am'],20./2]
    CS = ax.contour(xs, ys, A2, values_cont, colors='g', linewidths = 1)#, linewidths = 2
    cs_quiver = ax.contour(xs, ys, attr_sense_field/np.amax(attr_sense_field), [0.8], colors='y', linewidths = 1)

    # Extract contour line data
    quiver_xs=[]
    quiver_ys=[]
    quiver_us=[]
    quiver_vs=[]
    quiver_cols=[]
    for collection in cs_quiver.collections:
        print (collection.get_paths())
        for path in collection.get_paths():
            vertices = path.vertices
            print(path.vertices)
            print(path.vertices.size)
            x_coords = vertices[::16, 0]
            y_coords = vertices[::16, 1]
            
            # ~ x = xx_quiver[0,:]
            # ~ y = yy_quiver[:,1]
            x = xx[0,:]
            y = yy[:,1]
            print(x)
    
            # Find closest grid indices
            x_indices = np.round((x_coords - x[0]) / (x[1] - x[0])).astype(int)
            y_indices = np.round((y_coords - y[0]) / (y[1] - y[0])).astype(int)
    
            # Now you have the list of indices corresponding to the 2D line
            print(x_indices, y_indices)
            print (xx[y_indices,x_indices])
            quiver_xs.extend(xx[y_indices,x_indices])
            quiver_ys.extend(yy[y_indices,x_indices])
            quiver_us.extend(chemo_flux[1][y_indices,x_indices])
            quiver_vs.extend(chemo_flux[0][y_indices,x_indices])
            quiver_cols.extend(np.where(B_tot[y_indices,x_indices]>0 , E[y_indices,x_indices]/B_tot[y_indices,x_indices], 0.))

    
    

    colormap = cmr.get_sub_cmap('cmr.lilac_r', 0.2, 0.9)
    q = ax.quiver(quiver_xs, quiver_ys, quiver_us,quiver_vs,color=colormap(quiver_cols),units='inches',width=0.01, angles='xy', scale_units='xy', scale=0.000035*np.amax(B_tot), alpha=0.7)# ,units='xy'

    
    ax.set_aspect('equal', adjustable='box')


    if idx_lapse==0:
        scalebar = AnchoredSizeBar(ax.transData,
                               25, '25 mm', 'lower left', 
                               pad=0.2,
                               color='black',
                               frameon=False,
                               size_vertical=0.3, fontproperties=fontprops)
        ax.add_artist(scalebar)
        
                
        # Create a single large arrow as the color scale
        arrow_length = 60  # Length of the arrow for the color scale
        x_pos, y_pos = -50, 50  # Starting position of the arrow (axes coordinates)
        segments = 100
        arrow_segments = []
        
        for i in range(segments):
            segment_start = x_pos + (i / segments) * arrow_length
            segment_end = x_pos + ((i + 1) / segments) * arrow_length
            color = colormap(i / segments)
            arrow_segments.append(
                FancyArrow(
                    segment_start, y_pos, segment_end - segment_start, 0, 
                    width=1.8, color=color, length_includes_head=False, head_width=4.5 if i == segments - 1 else 0
                )
            )
        
        # Add colored arrow to the plot
        patch_collection = PatchCollection(arrow_segments, match_original=True)
        ax.add_collection(patch_collection)
        
        # Add labels for the quiver key
        ax.text(x_pos, y_pos - 8, f'{0:.1f}', ha='center', va='center', fontsize=17)
        ax.text(x_pos + arrow_length, y_pos - 8, f'{1:.1f}', ha='center', va='center', fontsize=17)
        ax.text(x_pos + arrow_length / 2, y_pos + 6, 'Infected fraction', ha='center', va='center', fontsize=17)
        
    ax.annotate(str(int(t))+" hr", xy=(35, 48), xycoords='data',xytext=(0, 5), textcoords='offset pixels', color='k', fontsize=17)             


    # ~ ax.legend(frameon=False, columnspacing=0.5, handletextpad=0.2,
      # ~ loc='upper right', bbox_to_anchor=(2.5, 1.2))
      
    ax.tick_params(left = False, right = False , labelleft = False ,
                labelbottom = False, bottom = False)
                
                
    return [heatmap_int,heatmap_V,heatmap_E,q]
    
fontprops = fm.FontProperties(size=17)

dir_fig='../../sims_for_fig5' # directory with input files
dir_io='{inp}/merge_T4_Vc1e7_60deg'.format(inp=dir_fig)  # directory with input files

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
t_int=[38]
t_quivers=[22,24]

figquivers, axs_figquivers = plt.subplots(nrows=1, ncols=len(t_quivers)+3, figsize=(11,4), width_ratios=[1, 1, 0.05,0.05,0.05], layout='constrained')# , layout='constrained',





for it, t in enumerate(times):
    
    if int(t) in t_int + t_ref + t_quivers:

        
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
        
        
        # PLOT INITIAL CONDITIONS
        
        in_x_virs= int((25000  -L/2.)/dx)*dx
        in_x_bacts= in_x_virs-(15000*np.sqrt(2)) # initial position of virs + 1cm shift, always microns

    
        in_y_bacts= 0# initial position of virs, always microns
        in_y_virs= 0# initial position of virs, always microns
        
        
        xx_rot_tmp_vir, yy_rot_tmp_vir = rotate((0,0),in_x_virs, in_y_virs, np.deg2rad(-45))
        xx_rot_tmp_bacts, yy_rot_tmp_bacts = rotate((0,0),in_x_bacts, in_y_bacts, np.deg2rad(-45))
        xx_rot_pivot, yy_rot_pivot = rotate((0,0),70000  -L/2.,0, np.deg2rad(-45))
    
        # ~ xx_rot, yy_rot = rotate((76000  -L/2.,0),xx_rot_tmp, np.abs(yy_rot_tmp), np.deg2rad(22.5))
        xx_rot_vir, yy_rot_vir = rotate((xx_rot_pivot,yy_rot_pivot),xx_rot_tmp_vir, yy_rot_tmp_vir, np.deg2rad(30))
        xx_rot_bacts, yy_rot_bacts = rotate((xx_rot_pivot,yy_rot_pivot),xx_rot_tmp_bacts, yy_rot_tmp_bacts, np.deg2rad(30))
        
        
              
        ax.plot(xx_rot_vir/1000, yy_rot_vir/1000, marker='x', markersize=8, color='k', linestyle='')
        ax.plot(xx_rot_bacts/1000, yy_rot_bacts/1000, marker='o',  markerfacecolor='none', markersize=8,   color='k', linestyle='')
        
        xx_rot_vir, yy_rot_vir = rotate((xx_rot_pivot,yy_rot_pivot),xx_rot_tmp_vir, yy_rot_tmp_vir, np.deg2rad(-30))
        xx_rot_bacts, yy_rot_bacts = rotate((xx_rot_pivot,yy_rot_pivot),xx_rot_tmp_bacts, yy_rot_tmp_bacts, np.deg2rad(-30))
        
        
              
        ax.plot(xx_rot_vir/1000, yy_rot_vir/1000, marker='x', markersize=8, color='k', linestyle='')
        ax.plot(xx_rot_bacts/1000, yy_rot_bacts/1000, marker='o',  markerfacecolor='none', markersize=8,   color='k', linestyle='')
        
  
        
        out_file='{out}/last_intensity.png'.format(out=dir_out_plots_tot)
        #    print out_file
        fig.savefig(out_file, dpi=300)
        fig.clf()
        plt.close('all')
        
    if int(t) in t_quivers:
        
        i_lapse= t_quivers.index(int(t)) 

           
        heatmap_int,heatmap_V,heatmap_E,q_fig2= plot_quiver(xx,yy,intensities,axs_figquivers[i_lapse],i_lapse)
            
            
        
 
 
 



cbar = figquivers.colorbar(heatmap_int, cax=axs_figquivers[-3])

cbar.ax.get_yaxis().labelpad = 20
cbar.ax.tick_params(axis ='y', labelsize=17)
cbar.ax.set_ylabel('pixel intensity', rotation=270, fontsize = 17)


cbar2 = figquivers.colorbar(heatmap_V,cax=axs_figquivers[-2])

cbar2.ax.get_yaxis().labelpad = 20
cbar2.ax.tick_params(axis ='y', labelsize=17)
cbar2.ax.set_ylabel(r"Phage density $(\rm{PFU} / mL)$", rotation=270, fontsize = 17)

cbar3 = figquivers.colorbar(heatmap_E, cax=axs_figquivers[-1])
cbar3.ax.get_yaxis().labelpad = 20
cbar3.ax.tick_params(axis ='y', labelsize=17)
cbar3.ax.set_ylabel(r"Infected hosts density $(\rm{CFU} / mL)$", rotation=270, fontsize = 17)


    

# ~ ax.set_aspect('equal', adjustable='box')


#### finish figure ####
labeldict = dict(labelstyle=r'{\sf \textbf{%s}}', fontsize='medium',
     xycoords=('axes fraction'), fontweight = 'bold')
#    mpsetup.label_axes([labeled_axes[0]], labels='A', xy=(-0.2, 0.95), **labeldict)
#mpsetup.label_axes([labeled_axes[1]], labels='B', xy=(-0.3, 0.95), **labeldict)
out_file='{out}/quiver_showcase.png'.format(out=dir_out_plots_tot)
#    print out_file
figquivers.savefig(out_file, dpi=300)
figquivers.clf()
plt.close('all')
    
 
    
 
    
dir_io='{inp}/merge_T4_Vc1e7_120deg'.format(inp=dir_fig)  # directory with input files

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
t_int=[54]
t_quivers=[20,22]

# ~ figquivers, axs_figquivers = plt.subplots(nrows=1, ncols=len(t_quivers), figsize=(11,4), gridspec_kw={'right':0.9, 'width_ratios':[1,1.5]})# , layout='constrained',
figquivers, axs_figquivers = plt.subplots(nrows=1, ncols=len(t_quivers)+3, figsize=(11,4), width_ratios=[1, 1, 0.05,0.05,0.05], layout='constrained')# , layout='constrained',
# ~ fig2, axs_fig2 = plt.subplots(nrows=3, ncols=len(times_quiverzoom)+2, figsize=(20,12),width_ratios=[1, 1,1,1, 0.05,0.05],height_ratios=[1, 1,1], layout='constrained')



    


for it, t in enumerate(times):
    
    if int(t) in t_int + t_ref + t_quivers:
        
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
        
        in_x_virs= int((25000  -L/2.)/dx)*dx
        in_x_bacts= in_x_virs-(15000*np.sqrt(2)) # initial position of virs + 1cm shift, always microns

    
        in_y_bacts= 0# initial position of virs, always microns
        in_y_virs= 0# initial position of virs, always microns
        
        
        xx_rot_tmp_vir, yy_rot_tmp_vir = rotate((0,0),in_x_virs, in_y_virs, np.deg2rad(-45))
        xx_rot_tmp_bacts, yy_rot_tmp_bacts = rotate((0,0),in_x_bacts, in_y_bacts, np.deg2rad(-45))
        
        xx_rot_pivot, yy_rot_pivot = rotate((0,0),70000  -L/2.,0, np.deg2rad(-45))
    
        # ~ xx_rot, yy_rot = rotate((76000  -L/2.,0),xx_rot_tmp, np.abs(yy_rot_tmp), np.deg2rad(22.5))
        xx_rot_vir, yy_rot_vir = rotate((xx_rot_pivot,yy_rot_pivot),xx_rot_tmp_vir, yy_rot_tmp_vir, np.deg2rad(60))
        xx_rot_bacts, yy_rot_bacts = rotate((xx_rot_pivot,yy_rot_pivot),xx_rot_tmp_bacts, yy_rot_tmp_bacts, np.deg2rad(60))
        
        
              
        ax.plot(xx_rot_vir/1000, yy_rot_vir/1000, marker='x', markersize=8, color='k', linestyle='')
        ax.plot(xx_rot_bacts/1000, yy_rot_bacts/1000, marker='o',  markerfacecolor='none', markersize=8,   color='k', linestyle='')
        
        xx_rot_vir, yy_rot_vir = rotate((xx_rot_pivot,yy_rot_pivot),xx_rot_tmp_vir, yy_rot_tmp_vir, np.deg2rad(-60))
        xx_rot_bacts, yy_rot_bacts = rotate((xx_rot_pivot,yy_rot_pivot),xx_rot_tmp_bacts, yy_rot_tmp_bacts, np.deg2rad(-60))
        
        
              
        ax.plot(xx_rot_vir/1000, yy_rot_vir/1000, marker='x', markersize=8, color='k', linestyle='')
        ax.plot(xx_rot_bacts/1000, yy_rot_bacts/1000, marker='o',  markerfacecolor='none', markersize=8,   color='k', linestyle='')
        
  
        
        out_file='{out}/last_intensity.png'.format(out=dir_out_plots_tot)
        #    print out_file
        fig.savefig(out_file, dpi=300)
        fig.clf()
        plt.close('all')
        
    if int(t) in t_quivers:
        
        i_lapse= t_quivers.index(int(t)) 
        
        heatmap_int,heatmap_V,heatmap_E,q_fig2= plot_quiver(xx,yy,intensities,axs_figquivers[i_lapse],i_lapse)
            
            

cbar = figquivers.colorbar(heatmap_int, cax=axs_figquivers[-3])

cbar.ax.get_yaxis().labelpad = 20
cbar.ax.tick_params(axis ='y', labelsize=17)
cbar.ax.set_ylabel('pixel intensity', rotation=270, fontsize = 17)


# ~ cax2 = divider.append_axes("right", size="5%", pad=0.9)
# ~ cbar2 = fig.colorbar(heatmap_V, cax=cax2)
cbar2 = figquivers.colorbar(heatmap_V,cax=axs_figquivers[-2])

cbar2.ax.get_yaxis().labelpad = 20
cbar2.ax.tick_params(axis ='y', labelsize=17)
cbar2.ax.set_ylabel(r"Phage density $(\rm{PFU} / mL)$", rotation=270, fontsize = 17)

# ~ cax3 = divider.append_axes("right", size="5%", pad=0.9)
# ~ cbar3 = fig.colorbar(heatmap_E, cax=cax3)
cbar3 = figquivers.colorbar(heatmap_E, cax=axs_figquivers[-1])
cbar3.ax.get_yaxis().labelpad = 20
cbar3.ax.tick_params(axis ='y', labelsize=17)
cbar3.ax.set_ylabel(r"Infected hosts density $(\rm{CFU} / mL)$", rotation=270, fontsize = 17)


    

# ~ ax.set_aspect('equal', adjustable='box')


#### finish figure ####
labeldict = dict(labelstyle=r'{\sf \textbf{%s}}', fontsize='medium',
     xycoords=('axes fraction'), fontweight = 'bold')
#    mpsetup.label_axes([labeled_axes[0]], labels='A', xy=(-0.2, 0.95), **labeldict)
#mpsetup.label_axes([labeled_axes[1]], labels='B', xy=(-0.3, 0.95), **labeldict)
out_file='{out}/quiver_showcase.png'.format(out=dir_out_plots_tot)
#    print out_file
figquivers.savefig(out_file, dpi=300)
figquivers.clf()
plt.close('all')
    
 
    

