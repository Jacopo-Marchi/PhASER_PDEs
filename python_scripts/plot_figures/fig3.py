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
    # PLOT INTENSITY SNAPSHOTS, ADAPTING LAYOUT TO THE SIZE OF THE INSET.
    v = ax.pcolormesh(xs, ys, zs, cmap='gray',vmin=0, vmax=255,shading='auto', linewidth=0,rasterized=True)
    ax.set_aspect('equal')
    
    

    if idx_lapse==0:
        scalebar = AnchoredSizeBar(ax.transData,
                               20, '2 cm', 'lower right', 
                               pad=0.2,
                               color='black',
                               frameon=False,
                               size_vertical=0.1, fontproperties=fontprops)
        ax.add_artist(scalebar)
      
    trans = transforms.blended_transform_factory(ax.transData, ax.get_xticklabels()[0].get_transform())
    if mfs_tit is not None:
        if '\n' not in mfs_tit:
            ax.annotate(mfs_tit, xy=(-90/rescdim, 85/rescdim), xycoords='data',xytext=(0, 5), textcoords='offset pixels', color='k', fontsize=16)
        else:
            ax.annotate(mfs_tit, xy=(-90/rescdim, 75/rescdim), xycoords='data',xytext=(0, 5), textcoords='offset pixels', color='k', fontsize=16)
        ax.annotate(str(int(t))+" hr", xy=(65/rescdim, 70/rescdim), xycoords='data',xytext=(0, 5), textcoords='offset pixels', color='k', fontsize=16)    #

    else:
        ax.annotate(str(int(t))+" hr", xy=(65/rescdim, 88/rescdim), xycoords='data',xytext=(0, 5), textcoords='offset pixels', color='k', fontsize=16)    #
    #
    ax.tick_params(left = False, right = False , labelleft = False ,
                labelbottom = False, bottom = False)
                
        
    if rightside is not None:
        width_inset=rightside-leftside
                    
        rect = patches.Rectangle((leftside, leftside+yoffset), width_inset, width_inset, linewidth=2, edgecolor='k', facecolor='none')

        # Add the patch to the Axes
        ax.add_patch(rect)
        
    return v

        
def plot_quiverzoom(xs,ys,zs,ax,idx_lapse,rightside=None,leftside=None,yoffset=None, rescdim=1):
    # PLOT INSETS WITH QUIVERS AND DETAILED INFORMATION ON THE SYSTEM STATES
    
    #PLOT intensities from all bacteria
    heatmap_int = ax.pcolormesh(xs, ys, zs, cmap='gray',vmin=0, vmax=255,shading='auto', linewidth=0,rasterized=True) # faster than pcolor
    
    cutvir=np.amax(V)/20. # introduce a cutoff to avoid doo much colors overlap
    V_cut= np.where(V>cutvir,V,np.nan)
    
    # plot phage
    heatmap_V = ax.pcolormesh(xs, ys, V_cut, cmap='Reds',vmin=1e6, vmax=1e10, alpha=0.5) # faster than pcolor, norm=colors.LogNorm()
  
    cutbact=np.amax(E)/20. # introduce a cutoff to avoid doo much colors overlap
    E_cut= np.where(E>cutbact,E,np.nan)
    
    # plot infected cells
    heatmap_E = ax.pcolormesh(xx, yy, E_cut, cmap='PuBu',vmin=1e6, vmax=5e8, alpha=0.5) # faster than pcolor, norm=colors.LogNorm()
        
    # plot two isolines for the attractant (aspartate) concentration at the lower sensing bound and half the initial concentration
    values_cont=[pars['am'],20./2]
    CS = ax.contour(xs, ys, A2, values_cont, colors='g', linewidths = 2)#, linewidths = 2
    
    # plot an isoline where the sensed gradient is at 80% its maximum, where the quivers are shown to represent the chemotactic flux vectors      
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
            x_coords = vertices[::8, 0] # sample a point every 8, otherwise quiver looks too dense 
            y_coords = vertices[::8, 1]
            
            x = xx[0,:]
            y = yy[:,1]
            print(x)
    
            # Find closest grid indices on the isoline
            x_indices = np.round((x_coords - x[0]) / (x[1] - x[0])).astype(int)
            y_indices = np.round((y_coords - y[0]) / (y[1] - y[0])).astype(int)
            # Now you have the list of indices corresponding to the 2D line
            
            print(x_indices, y_indices)
            print (xx[y_indices,x_indices])
            # store sample coordinates and corresponding flux vectors to plot later
            quiver_xs.extend(xx[y_indices,x_indices])
            quiver_ys.extend(yy[y_indices,x_indices])
            quiver_us.extend(chemo_flux[1][y_indices,x_indices])
            quiver_vs.extend(chemo_flux[0][y_indices,x_indices])
            quiver_cols.extend(np.where(B_tot[y_indices,x_indices]>0 , E[y_indices,x_indices]/B_tot[y_indices,x_indices], 0.))

    
    

    colormap = cmr.get_sub_cmap('cmr.lilac_r', 0.2, 0.9)
    q = ax.quiver(quiver_xs, quiver_ys, quiver_us,quiver_vs,color=colormap(quiver_cols),units='inches',width=0.022, angles='xy', scale_units='xy', scale=0.00007*np.amax(B_tot))# ,units='xy'

    
    ax.set_aspect('equal', adjustable='box')

    if idx_lapse==0: # plot scale bar
        scalebar = AnchoredSizeBar(ax.transData,
                               5, '5 mm', 'upper left', 
                               pad=0.2,
                               color='black',
                               frameon=False,
                               size_vertical=0.1, fontproperties=fontprops)
        ax.add_artist(scalebar)
        
                
        # Create a single large arrow as the color scale
        arrow_length = 15/rescdim  # Length of the arrow for the color scale
        x_pos, y_pos = leftside+8/rescdim, rightside+yoffset-4/rescdim  # Starting position of the arrow (axes coordinates)
        segments = 100
        arrow_segments = []
        
        for i in range(segments):
            segment_start = x_pos + (i / segments) * arrow_length
            segment_end = x_pos + ((i + 1) / segments) * arrow_length
            color = colormap(i / segments)
            arrow_segments.append(
                FancyArrow(
                    segment_start, y_pos, segment_end - segment_start, 0, 
                    width=0.5, color=color, length_includes_head=False, head_width=1.3 if i == segments - 1 else 0
                )
            )
        
        # Add colored arrow to the plot
        patch_collection = PatchCollection(arrow_segments, match_original=True)
        ax.add_collection(patch_collection)
        
        # Add labels for the quiver key
        ax.text(x_pos, y_pos - 1.5, f'{0:.1f}', ha='center', va='center', fontsize=14)
        ax.text(x_pos + arrow_length, y_pos - 1.5, f'{1:.1f}', ha='center', va='center', fontsize=14)
        ax.text(x_pos + arrow_length / 2, y_pos + 1.2, 'Infected fraction', ha='center', va='center', fontsize=16)

      
    ax.tick_params(left = False, right = False , labelleft = False ,
                labelbottom = False, bottom = False)
                
    ax.annotate(str(int(t))+" hr", xy=(rightside-4/rescdim, rightside+yoffset-1.5/rescdim), xycoords='data',xytext=(0, 5), textcoords='offset pixels', color='k', fontsize=16)             
                
    ax.set_xlim(xmin=leftside, xmax=rightside)
    ax.set_ylim(ymin=leftside+yoffset, ymax=rightside+yoffset)
    
    return [heatmap_int,heatmap_V,heatmap_E,q]

    
    
fontprops = fm.FontProperties(size=16)

dir_fig='../../sims_for_fig3' # directory with input files
dir_io='{inp}/T7_fin_Vcinf'.format(inp=dir_fig)  # directory with input files

dir_in_tot='{inp}/data'.format(inp=dir_io) # 
dir_in_frames='{inp}/frames'.format(inp=dir_in_tot) # 
dir_out_plots_tot='{inp}/plots'.format(inp=dir_io) # 

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



# load calibration from CFu to intensity
df = pd.read_excel('{inp}/../../growth_data/Scanner_OD.xlsx'.format(inp=dir_io), skiprows=1, usecols='A,B',  header = None,  engine='openpyxl', sheet_name="CFU Calibration")


pxtoCFUs=df.to_numpy(copy=True)[:8,:]

pxtoCFUs=pxtoCFUs[np.argsort(pxtoCFUs[:,1]),:]

print(pxtoCFUs)
print(pxtoCFUs[0,:])

f = make_interp_spline(pxtoCFUs[:,1], pxtoCFUs[:,0], k=3)


       

# times of snapshots and zoom borders, in data coordinates
times_lapses=[24, 38, 40, 56]
times_quiverzoom=times_lapses
leftside_zoom=[-35, -7, -3, 60]
rightside_zoom=[-10, 18, 22, 85]
yoffset_zoom=[-23, 0, 0, 0]



fig_timelapse, axs_timelapse = plt.subplots(nrows=1, ncols=len(times_lapses), figsize=(18,4))
fig_quiverzoom, axs_quiverzoom = plt.subplots(nrows=1, ncols=len(times_quiverzoom), figsize=(21,4),width_ratios=[1, 1,1,1.5])
fig2, axs_fig2 = plt.subplots(nrows=3, ncols=len(times_quiverzoom)+2, figsize=(20,12),width_ratios=[1, 1,1,1, 0.05,0.05],height_ratios=[1, 1,1], layout='constrained')

    
labeled_axes = []
labeled_axes.append(axs_fig2[0,0])
labeled_axes.append(axs_fig2[2,0])
labeled_axes.append(axs_fig2[2,1])
labeled_axes.append(axs_fig2[2,2])
    

# load data for PANEL A
savedat=True
for file_in_frames in glob.glob("{inp}/experiment_state_time_*".format(inp=dir_in_frames)):
    filename, ext = os.path.splitext(file_in_frames)
    savedat=(ext=='.dat')
    print(ext, savedat)
    time=os.path.basename(filename).split('_')[-1]
    print(time)
    # ~ time=int(float(time))
    t=float(time)

 
    print(t)

    if int(t) in times_lapses:
        
        if savedat:
            data_space = np.loadtxt(file_in_frames)
        else:
            with np.load(file_in_frames) as data:
                 data_space = data['arr_0']
     
    
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
    
        B_tot = S+E
        
            
        derivs_plots = derivatives_euler(pars, 'derivatives')
        
        newstate = derivs_plots(R,S,V,E, A,A2)  # compute derivatives the same way it's done during PDE integration'
        
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

        i_lapse= times_lapses.index(int(t)) 
        
        v= plot_timelapse_intensity(xx,yy,intensities,axs_timelapse[i_lapse],i_lapse,rightside_zoom[i_lapse], leftside_zoom[i_lapse],yoffset_zoom[i_lapse]) #,rescdim=10/6
        
        heatmap_int,heatmap_V,heatmap_E,q= plot_quiverzoom(xx,yy,intensities,axs_quiverzoom[i_lapse],i_lapse,rightside_zoom[i_lapse], leftside_zoom[i_lapse],yoffset_zoom[i_lapse])
        
        v_fig2= plot_timelapse_intensity(xx,yy,intensities,axs_fig2[0,i_lapse],i_lapse,rightside_zoom[i_lapse], leftside_zoom[i_lapse],yoffset_zoom[i_lapse]) # plots first row of Panel A
        
        heatmap_int_fig2,heatmap_V_fig2,heatmap_E_fig2,q_fig2= plot_quiverzoom(xx,yy,intensities,axs_fig2[1,i_lapse],i_lapse,rightside_zoom[i_lapse], leftside_zoom[i_lapse],yoffset_zoom[i_lapse]) # plots second row of Panel A with quiver insets
        
        
                
            


divider = make_axes_locatable(axs_timelapse[len(axs_timelapse)-1])
cax1 = divider.append_axes("right", size="5%", pad=0.2)
cbar = fig_timelapse.colorbar(v, cax=cax1)
cbar.ax.get_yaxis().labelpad = 20
cbar.ax.tick_params(axis ='y', labelsize=14)
cbar.ax.set_ylabel('pixel intensity', rotation=270, fontsize = 16)

fig_timelapse.tight_layout()
out_file='{out}/wave_timelapse_sep_panels.pdf'.format(out=dir_out_plots_tot, time=t)
fig_timelapse.savefig(out_file)



divider = make_axes_locatable(axs_quiverzoom[len(axs_quiverzoom)-1])
cax1 = divider.append_axes("right", size="5%", pad=0.2)
cbar = fig_quiverzoom.colorbar(heatmap_int, cax=cax1)
cbar.ax.get_yaxis().labelpad = 20
cbar.ax.tick_params(axis ='y', labelsize=14)
cbar.ax.set_ylabel('pixel intensity', rotation=270, fontsize = 16)

cax2 = divider.append_axes("right", size="5%", pad=0.9)
cbar2 = fig_quiverzoom.colorbar(heatmap_V, cax=cax2)
cbar2.ax.get_yaxis().labelpad = 20
cbar2.ax.tick_params(axis ='y', labelsize=14)
cbar2.ax.set_ylabel(r"Phage density $(\rm{PFU} / mL)$", rotation=270, fontsize = 16)

cax3 = divider.append_axes("right", size="5%", pad=0.9)
cbar3 = fig_quiverzoom.colorbar(heatmap_E, cax=cax3)
cbar3.ax.get_yaxis().labelpad = 35
cbar3.ax.tick_params(axis ='y', labelsize=14)
cbar3.ax.set_ylabel("Infected bacteria \n"+r"density $(\rm{CFU} / mL)$", rotation=270, fontsize = 16)


fig_quiverzoom.tight_layout()
out_file='{out}/quiverzoom_sep_panels.png'.format(out=dir_out_plots_tot)
fig_quiverzoom.savefig(out_file, dpi=300)


# plot color bars for panel A
divider = make_axes_locatable(axs_fig2[0,len(axs_timelapse)-1])
cbar = fig2.colorbar(v_fig2, cax=axs_fig2[0,-2])
cbar.ax.get_yaxis().labelpad = 20
cbar.ax.tick_params(axis ='y', labelsize=14)
cbar.ax.set_ylabel('pixel intensity', rotation=270, fontsize = 16)


cbar2 = fig2.colorbar(heatmap_V_fig2, cax=axs_fig2[1,-2])
cbar2.ax.get_yaxis().labelpad = 20
cbar2.ax.tick_params(axis ='y', labelsize=14)
cbar2.ax.set_ylabel(r"Phage density $(\rm{PFU} / mL)$", rotation=270, fontsize = 16)

cbar3 = fig2.colorbar(heatmap_E_fig2, cax=axs_fig2[1,-1], pad=0.04)
cbar3.ax.get_yaxis().labelpad = 20
cbar3.ax.tick_params(axis ='y', labelsize=14)
cbar3.ax.set_ylabel(r"Infected hosts density $(\rm{CFU} / mL)$", rotation=270, fontsize = 16)

fig2.text(0.005,0.5,'Phage T7, single attractant', rotation='vertical', fontsize=16)

axs_fig2[0,-1].set_visible(False)
axs_fig2[2,-2].set_visible(False)
axs_fig2[2,-1].set_visible(False)


####PANEL B, two attractants

dir_io='{inp}/T7_fin_Vcinf_F8'.format(inp=dir_fig)  # directory with input files

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

t_ref=[50] # time of plotted snapshot
savedat=True
for file_in_frames in glob.glob("{inp}/experiment_state_time_*".format(inp=dir_in_frames)):
    filename, ext = os.path.splitext(file_in_frames)
    savedat=(ext=='.dat')
    print(ext, savedat)
    time=os.path.basename(filename).split('_')[-1]
    print(time)
    # ~ time=int(float(time))
    t=float(time)

 
    print(t)

    if int(t) in t_ref:
        
        if savedat:
            data_space = np.loadtxt(file_in_frames)
        else:
            with np.load(file_in_frames) as data:
                 data_space = data['arr_0']
     
    
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
        
    
        xnew = B_tot
        intensities = f(xnew)
        intensities = np.array(intensities).astype(float)  
        
        print(intensities)
        print(intensities.dtype)
        
        
        xx = np.reshape(x, (Nrows, Ncols))/1000.
        print (xx.shape)
        yy = np.reshape(y, xx.shape)/1000.
        R = np.reshape(R, xx.shape)
        S = np.reshape(S, xx.shape)
        E = np.reshape(E, xx.shape)
        V = np.reshape(V, xx.shape)
        A = np.reshape(A, xx.shape)
        A2 = np.reshape(A2, xx.shape)
        B_tot = np.reshape(B_tot, xx.shape)
    
        intensities = np.reshape(intensities, xx.shape)
    
        B_tot = S+E
        
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
    
            
        
        i_lapse= t_ref.index(int(t)) 

        v_fig2= plot_timelapse_intensity(xx,yy,intensities,axs_fig2[2,i_lapse],i_lapse, mfs_tit='Phage T7, double attractant')
        
        
        
                
    

####PANEL C, T4

dir_io='{inp}/T4_fin_Vc1e7_F8'.format(inp=dir_fig)  # directory with input files

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

t_ref=[50]
savedat=True
for file_in_frames in glob.glob("{inp}/experiment_state_time_*".format(inp=dir_in_frames)):
    filename, ext = os.path.splitext(file_in_frames)
    savedat=(ext=='.dat')
    print(ext, savedat)
    time=os.path.basename(filename).split('_')[-1]
    print(time)
    # ~ time=int(float(time))
    t=float(time)

 
    print(t)

    if int(t) in t_ref:
        
        if savedat:
            data_space = np.loadtxt(file_in_frames)
        else:
            with np.load(file_in_frames) as data:
                 data_space = data['arr_0']
     
    
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
        
        xnew = B_tot
        intensities = f(xnew)
        intensities = np.array(intensities).astype(float)  
        
        print(intensities)
        print(intensities.dtype)
        
        
        xx = np.reshape(x, (Nrows, Ncols))/1000.
        print (xx.shape)
        yy = np.reshape(y, xx.shape)/1000.
        R = np.reshape(R, xx.shape)
        S = np.reshape(S, xx.shape)
        E = np.reshape(E, xx.shape)
        V = np.reshape(V, xx.shape)
        A = np.reshape(A, xx.shape)
        A2 = np.reshape(A2, xx.shape)
        B_tot = np.reshape(B_tot, xx.shape)
    
        intensities = np.reshape(intensities, xx.shape)
    
        B_tot = S+E
        
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
        
        i_lapse= t_ref.index(int(t)) 

        v_fig2= plot_timelapse_intensity(xx,yy,intensities,axs_fig2[2,i_lapse+1],i_lapse, mfs_tit='Phage T4, double attractant')#,rescdim=10/6
        
    
    

####PANEL D, T7, with slower infected cells, composed of two plots

dir_in_mfs='{inp}/modulate_chemoeff_t7/'.format(inp=dir_fig)  # directory with input files

mfs=[0.7,0.75] # motility difference factors. The simulated data need to be in the target folder
mfs_str=['0d7','0d75']
mfs_tits=['Phage T7, 70 % infected \n chemosensing','Phage T7, 75 % infected \n chemosensing']

for i_mf,mf in enumerate(mfs):
    dir_io ="{inp}/t7_FIN_mf_{mf}".format(inp=dir_in_mfs,mf=mfs_str[i_mf])
    
    
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

    t_ref=[40]
    savedat=True
    for file_in_frames in glob.glob("{inp}/experiment_state_time_*".format(inp=dir_in_frames)):
        filename, ext = os.path.splitext(file_in_frames)
        savedat=(ext=='.dat')
        print(ext, savedat)
        time=os.path.basename(filename).split('_')[-1]
        print(time)
        # ~ time=int(float(time))
        t=float(time)
    
     
        print(t)
    
        if int(t) in t_ref:
            
            if savedat:
                data_space = np.loadtxt(file_in_frames)
            else:
                with np.load(file_in_frames) as data:
                     data_space = data['arr_0']
         
        
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
            
            xnew = B_tot
            intensities = f(xnew)
            intensities = np.array(intensities).astype(float)  
            
            print(intensities)
            print(intensities.dtype)
            
            
            xx = np.reshape(x, (Nrows, Ncols))/1000.
            print (xx.shape)
            yy = np.reshape(y, xx.shape)/1000.
            R = np.reshape(R, xx.shape)
            S = np.reshape(S, xx.shape)
            E = np.reshape(E, xx.shape)
            V = np.reshape(V, xx.shape)
            A = np.reshape(A, xx.shape)
            A2 = np.reshape(A2, xx.shape)
            B_tot = np.reshape(B_tot, xx.shape)
        
            intensities = np.reshape(intensities, xx.shape)
        
            B_tot = S+E
                
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
        
                
            
            i_lapse= t_ref.index(int(t)) 
    
            v_fig2= plot_timelapse_intensity(xx,yy,intensities,axs_fig2[2,i_lapse+i_mf+2],i_lapse+i_mf, mfs_tit=mfs_tits[i_mf])
            
            
                    
        
    
    


#### finish figure ####
labeldict = dict(labelstyle=r'%s', fontsize=26,
     xycoords=('axes fraction'), fontweight = 'bold')
mpsetup.label_axes([labeled_axes[0]], labels='A', xy=(-0.1, 0.95), **labeldict)
mpsetup.label_axes([labeled_axes[1]], labels='B', xy=(-0.1, 0.95), **labeldict)
mpsetup.label_axes([labeled_axes[2]], labels='C', xy=(-0.1, 0.95), **labeldict)
mpsetup.label_axes([labeled_axes[3]], labels='D', xy=(-0.1, 0.95), **labeldict)
    

out_file='{out}/Fig3.png'.format(out=dir_fig)

fig2.savefig(out_file, dpi=300)


