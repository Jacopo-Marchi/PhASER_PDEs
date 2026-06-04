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
                               25, '25 mm', 'lower right', 
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

        
    
fontprops = fm.FontProperties(size=16)

dir_fig='../../sims_for_fig4_new' # directory with input files



df = pd.read_excel('{inp}/../growth_data/Scanner_OD.xlsx'.format(inp=dir_fig), skiprows=1, usecols='A,B',  header = None,  engine='openpyxl', sheet_name="CFU Calibration")


pxtoCFUs=df.to_numpy(copy=True)[:8,:]

pxtoCFUs=pxtoCFUs[np.argsort(pxtoCFUs[:,1]),:]

print(pxtoCFUs)
print(pxtoCFUs[0,:])

# ~ f = interpolate.interp1d(pxtoCFUs[:,1], pxtoCFUs[:,0])
f = make_interp_spline(pxtoCFUs[:,1], pxtoCFUs[:,0], k=3)

############ PLOT FRAMES



dirstrings=['{inp}/F8_t4_Vc1e7_mf0d87_P05e5','{inp}/t4_F8_mf1_P05e5','{inp}/t7_F8_mf0d87_P05e5','{inp}/F8_t7_Vcinf_P05e5']
t_ints=[26,44,44, 42]

figallee, axs_allee = plt.subplots(nrows=2, ncols=(len(dirstrings)+1)//2, figsize=(8,8),  layout='constrained')# , layout='constrained',

# ~ chemoeffs = [
    # ~ "87%",
    # ~ "100%"
# ~ ]

chemoeffs = [
    "Infected chemosensing 87%",
    "Infected chemosensing 100%"
]
    
labeled_axes = []
labeled_axes.append(axs_allee[0,0])
# ~ labeled_axes.append(axs_fig2[1,0])
labeled_axes.append(axs_allee[0,1])
labeled_axes.append(axs_allee[1,0])
labeled_axes.append(axs_allee[1,1])


# ~ fig_ins = plt.figure(figsize=(thisfigsize[0]*4./3, thisfigsize[0]))
# ~ grid_ins = gridspec.GridSpec(2, 2, left=0.05, right=0.93, top=0.95, bottom=0.07,
     # ~ wspace=0.35, hspace=0.4)
# ~ ax_ins = plt.Subplot(fig_ins, grid_ins[1, 0])#, rasterized=True
# ~ fig_ins.add_subplot(ax_ins)
    
for idir, dirstr in enumerate(dirstrings):


    
    
    dir_io=dirstr.format(inp=dir_fig)  # directory with input files
    
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
    
    
    bactprec=None    
    intprec=None    
    wavepos=None    
    waveposprec=None    
    frontpos=None    
    frontposprec=None    
    wavesum=0    
    
    
    t_estPFUs=[8]
    t_int=t_ints[idir]
    
    PFU_est=0
    
    
    savedat=True
    for file_in_frames in glob.glob("{inp}/experiment_state_time_*".format(inp=dir_in_frames)):
        filename, ext = os.path.splitext(file_in_frames)
        savedat=(ext=='.dat')
        print(ext, savedat)
        time=os.path.basename(filename).split('_')[-1]
        print(time)
        # ~ time=int(float(time))
        t=float(time)
    
      
        
        if int(t) in  t_estPFUs:
    
            
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
                
            
            # ~ PFU_est= np.sum(V[V/discr_thr>1]/discr_thr)
            PFU_est= np.amax(V/discr_thr)

        
    for file_in_frames in glob.glob("{inp}/experiment_state_time_*".format(inp=dir_in_frames)):
        filename, ext = os.path.splitext(file_in_frames)
        savedat=(ext=='.dat')
        print(ext, savedat)
        time=os.path.basename(filename).split('_')[-1]
        print(time)
        # ~ time=int(float(time))
        t=float(time)
    
      
        
        if int(t) == t_int :
    
            
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
                
                 
            row_plot= int(idir/2) 
            col_plot= np.mod(idir,2)
                
            v= plot_timelapse_intensity(xx,yy,intensities,axs_allee[row_plot,col_plot],0) #,rescdim=10/6
            
            
            axs_allee[row_plot, col_plot].plot((21000 -L/2.)/1000, (21000 -L/2.)/1000, marker='x', markersize=6, color='k', linestyle='')
            axs_allee[row_plot, col_plot].plot((6000  -L/2.)/1000, (6000  -L/2.)/1000, marker='o',  markerfacecolor='none', markersize=6,   color='k', linestyle='')
            
            if not row_plot:
                axs_allee[row_plot, col_plot].set_title(chemoeffs[col_plot], fontsize=16) # , y=-0.1

            
            if L > 120000:
                axs_allee[row_plot, col_plot].set_xlim(xmin=-100, xmax= 20)
                axs_allee[row_plot, col_plot].set_ylim(ymin=-100, ymax= 20)

            
                     
row_labels = [
    "T4",
    "T7"
]

for row_idx, label in enumerate(row_labels):
    axs_allee[row_idx, 0].annotate(
        label,
        xy=(0, 0.5),
        xycoords='axes fraction',
        xytext=(-12, -2),
        textcoords='offset points',
        ha='right',
        va='center',
        fontsize=16
    )        

# ~ axs_allee[0, 0].annotate(
        # ~ "Infected chemosensing:",
        # ~ xy=(0.4, 1.05),
        # ~ xycoords='axes fraction',
        # ~ xytext=(12, -1),
        # ~ textcoords='offset points',
        # ~ ha='right',
        # ~ va='center',
        # ~ fontsize=16
    # ~ )        
# ~ bbox = axs_allee[1,1].get_window_extent().transformed(figallee.dpi_scale_trans.inverted())
# ~ width, height = bbox.width, bbox.height

# ~ width *= figallee.dpi
# ~ height *= figallee.dpi

# ~ axs_allee[0,0].set_visible(False)

# ~ axs_allee[1, 0].set_position([-0.14, 0.26, 0.48, 0.48])

# ~ ax.set_aspect('equal', adjustable='box')

#### finish figure ####
labeldict = dict(labelstyle=r'%s', fontsize=26,
     xycoords=('axes fraction'), fontweight = 'bold')
mpsetup.label_axes([labeled_axes[0]], labels='A', xy=(-0.1, 0.93), **labeldict)
mpsetup.label_axes([labeled_axes[1]], labels='B', xy=(-0.1, 0.93), **labeldict)
mpsetup.label_axes([labeled_axes[2]], labels='C', xy=(-0.1, 0.93), **labeldict)
mpsetup.label_axes([labeled_axes[3]], labels='D', xy=(-0.1, 0.93), **labeldict)
out_file='{out}/SIfig_allee_diffchemo.png'.format(out=dir_fig)
#    print out_file
figallee.savefig(out_file, dpi=300)
figallee.clf()
plt.close('all')
    
 
    
 
