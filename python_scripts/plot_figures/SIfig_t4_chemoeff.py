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
                               20, '2 cm', 'lower right', 
                               pad=0.2,
                               color='black',
                               frameon=False,
                               size_vertical=0.1, fontproperties=fontprops)
        ax.add_artist(scalebar)
    # ~ ax.legend(frameon=False, columnspacing=0.5, handletextpad=0.2,
      # ~ loc='upper right', bbox_to_anchor=(2.5, 1.2))
      
    trans = transforms.blended_transform_factory(ax.transData, ax.get_xticklabels()[0].get_transform())
    if mfs_tit is not None:
        if '\n' not in mfs_tit:
            ax.annotate(mfs_tit, xy=(-90/rescdim, 85/rescdim), xycoords='data',xytext=(0, 5), textcoords='offset pixels', color='k', fontsize=16)
        else:
            ax.annotate(mfs_tit, xy=(-90/rescdim, 75/rescdim), xycoords='data',xytext=(0, 5), textcoords='offset pixels', color='k', fontsize=16)
        # ~ ax.annotate(str(int(t))+" hr", xy=(65/rescdim, 70/rescdim), xycoords='data',xytext=(0, 5), textcoords='offset pixels', color='k', fontsize=16)    #
        ax.annotate(str(int(t))+" hr", xy=(-90/rescdim, 55/rescdim), xycoords='data',xytext=(0, 5), textcoords='offset pixels', color='k', fontsize=16)    #

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

    
    step=6
    # ~ step=1
    # ~ xx_quiver=xx[::step,::step]
    # ~ yy_quiver=yy[::step,::step]
    # ~ chemo_flux_x_quiver=chemo_flux[1][::step,::step]
    # ~ chemo_flux_y_quiver=chemo_flux[0][::step,::step]
    # ~ grad_sensing_x_quiver=grad_sensing[1][::step,::step]
    # ~ grad_sensing_y_quiver=grad_sensing[0][::step,::step]

    # ~ bool_array= (attr_sense_field > 0.6*np.amax(attr_sense_field)) # divides sensing field below and above 70%
    # ~ bool_array_step = (attr_sense_field[::step,::step] > 0.7*np.amax(attr_sense_field[::step,::step])) # divides sensing field below and above 70%
    # ~ mask_quiver=  (np.diff(bool_array, axis=0, prepend=bool_array[0,0])) | (np.diff(bool_array, prepend=bool_array[0,0]) ) # check where difference is above zero (crossing 70%) in either direction
    # ~ mask_quiver_step=  (np.diff(bool_array_step, axis=0, prepend=bool_array_step[0,0])) | (np.diff(bool_array_step, prepend=bool_array_step[0,0]) ) # check where difference is above zero (crossing 70%) in either direction

    # ~ mask_quiver = mask_quiver & (B_tot> np.amax(B_tot)/20) # I also want bacteria to be at least half max 
    
    
    # ~ thr_attr=np.percentile(attr_sense_field[::step,::step][mask_quiver_step], 50)
    # ~ mask_quiver_step= mask_quiver_step & (attr_sense_field[::step,::step]>thr_attr)


    # ~ mask_quiver=  mask_quiver[::step,::step]
    

    # ~ heatmap = ax.pcolormesh(xx, yy, intensities, cmap='gray',vmin=maxbacts_convert, vmax=no_bactcontrol) # faster than pcolor
    heatmap_int = ax.pcolormesh(xs, ys, zs, cmap='gray',vmin=0, vmax=255,shading='auto', linewidth=0,rasterized=True) # faster than pcolor
    
    cutvir=np.amax(V)/20.
    V_cut= np.where(V>cutvir,V,np.nan)
    

    heatmap_V = ax.pcolormesh(xs, ys, V_cut, cmap='Reds',vmin=1e6, vmax=1e10, alpha=0.5) # faster than pcolor, norm=colors.LogNorm()
  
    # ~ cutbact=min(discr_thr, np.amax(E)/20.)
    cutbact=np.amax(E)/20.
    E_cut= np.where(E>cutbact,E,np.nan)

    heatmap_E = ax.pcolormesh(xx, yy, E_cut, cmap='PuBu',vmin=1e6, vmax=5e8, alpha=0.5) # faster than pcolor, norm=colors.LogNorm()
        
    values_cont=[pars['am'],20./2]
    CS = ax.contour(xs, ys, A2, values_cont, colors='g', linewidths = 2)#, linewidths = 2
            
    # ~ grad_modulo= np.sqrt(chemo_flux[1]**2. + chemo_flux[0]**2.)
    # ~ grad_modulo= np.where(E>cutbact,grad_modulo, 0.)

    
    
    # ~ col=(E[::step,::step]/B_tot[::step,::step])
    # ~ col=np.where(B_tot[::step,::step]>0 , col, 0.)
            
    
    # ~ cs_quiver = axs_quiverzoom[i_lapse].contour(xx_quiver, yy_quiver, attr_sense_field[::step,::step]/np.amax(attr_sense_field[::step,::step]), [0.7], colors='y', linewidths = 1)
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
            x_coords = vertices[::8, 0]
            y_coords = vertices[::8, 1]
            
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
            # ~ idx=np.vstack([x_indices,y_indices])
            # ~ print(idx)
            # ~ list_idxs=tuple(idx.T)
            # ~ list_idxs=list(zip(x_indices,y_indices))
            # ~ print(list_idxs)
            # ~ print (xx_quiver[y_indices,x_indices])
            # ~ quiver_xs.extend(xx_quiver[y_indices,x_indices])
            # ~ quiver_ys.extend(yy_quiver[y_indices,x_indices])
            # ~ quiver_us.extend(chemo_flux_x_quiver[y_indices,x_indices])
            # ~ quiver_vs.extend(chemo_flux_y_quiver[y_indices,x_indices])
            # ~ quiver_cols.extend(col[y_indices,x_indices])
            print (xx[y_indices,x_indices])
            quiver_xs.extend(xx[y_indices,x_indices])
            quiver_ys.extend(yy[y_indices,x_indices])
            quiver_us.extend(chemo_flux[1][y_indices,x_indices])
            quiver_vs.extend(chemo_flux[0][y_indices,x_indices])
            quiver_cols.extend(np.where(B_tot[y_indices,x_indices]>0 , E[y_indices,x_indices]/B_tot[y_indices,x_indices], 0.))

    # ~ print (quiver_xs)        
    # ~ print (quiver_us) 
    # ~ print (quiver_cols) 
           
    # ~ col=col[mask_quiver_step]        
    # ~ print (col)
    
    # ~ norm = colors.Normalize()
    # ~ norm.autoscale(quiver_cols)
    # we need to normalize our colors array to match it colormap domain
    # which is [0, 1]
    
    

    colormap = cmr.get_sub_cmap('cmr.lilac_r', 0.2, 0.9)
    # ~ q = ax.quiver(quiver_xs, quiver_ys, quiver_us,quiver_vs,color=colormap(norm(quiver_cols)),units='inches',width=0.02, angles='xy', scale_units='xy', scale=0.0001*np.amax(B_tot))# ,units='xy'
    q = ax.quiver(quiver_xs, quiver_ys, quiver_us,quiver_vs,color=colormap(quiver_cols),units='inches',width=0.022, angles='xy', scale_units='xy', scale=0.00007*np.amax(B_tot))# ,units='xy'

    
    ax.set_aspect('equal', adjustable='box')


    if idx_lapse==0:
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

    # ~ ax.legend(frameon=False, columnspacing=0.5, handletextpad=0.2,
      # ~ loc='upper right', bbox_to_anchor=(2.5, 1.2))
      
    ax.tick_params(left = False, right = False , labelleft = False ,
                labelbottom = False, bottom = False)
                
    ax.annotate(str(int(t))+" hr", xy=(rightside-4/rescdim, rightside+yoffset-1.5/rescdim), xycoords='data',xytext=(0, 5), textcoords='offset pixels', color='k', fontsize=16)             
                
    ax.set_xlim(xmin=leftside, xmax=rightside)
    ax.set_ylim(ymin=leftside+yoffset, ymax=rightside+yoffset)
    
    return [heatmap_int,heatmap_V,heatmap_E,q]

    
    
fontprops = fm.FontProperties(size=16)

dir_fig='../../sims_for_fig1_quiversinglefront' # directory with input files
dir_io='{inp}/T7_fin_Vcinf'.format(inp=dir_fig)  # directory with input files

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

       


fig2, axs_fig2 = plt.subplots(nrows=1, ncols=2, figsize=(10,4),layout='constrained')

# ~ fig_ins = plt.figure(figsize=(thisfigsize[0]*4./3, thisfigsize[0]))
# ~ grid_ins = gridspec.GridSpec(2, 2, left=0.05, right=0.93, top=0.95, bottom=0.07,
     # ~ wspace=0.35, hspace=0.4)
# ~ ax_ins = plt.Subplot(fig_ins, grid_ins[1, 0])#, rasterized=True
# ~ fig_ins.add_subplot(ax_ins)
    
labeled_axes = []
# ~ labeled_axes.append(axs_fig2[0,0])
# ~ labeled_axes.append(axs_fig2[2,0])
# ~ labeled_axes.append(axs_fig2[2,1])
# ~ labeled_axes.append(axs_fig2[2,2])
    


    

####PANEL E

dir_in_mfs='{inp}/modulate_chemoeff_t4/'.format(inp=dir_fig)  # directory with input files

mfs=[0.85,0.87]
mfs_str=['0d85','0d87']
mfs_tits=['Phage T4, 85 % infected \n chemosensing','Phage T4, 87 % infected \n chemosensing']

for i_mf,mf in enumerate(mfs):
# ~ for dir_io in glob.glob("{inp}/t7_FIN_mf_*".format(inp=dir_in_mfs)):
    dir_io ="{inp}/t4_VC1e7_FIN_huge_mf_{mf}".format(inp=dir_in_mfs,mf=mfs_str[i_mf])
    
    
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

    t_ref=[54]
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
        
                
            
            i_lapse= t_ref.index(int(t)) 
    
            v_fig2= plot_timelapse_intensity(xx,yy,intensities,axs_fig2[i_mf],0, mfs_tit=mfs_tits[i_mf])
            
            # ~ heatmap_int_fig2,heatmap_V_fig2,heatmap_E_fig2,q_fig2= plot_quiverzoom(xx,yy,intensities,axs_fig2[1,i_lapse],i_lapse,rightside_zoom[i_lapse], leftside_zoom[i_lapse],yoffset_zoom[i_lapse])
            
                    
        
    
    

# ~ fig2.tight_layout()
# ~ out_file='{out}/quiverzoom_sep_panels.pdf'.format(out=dir_out_plots_tot, time=t)
out_file='{out}/SIfig_t4_chemoeff.png'.format(out=dir_fig)

fig2.savefig(out_file, dpi=300)

#    print out_file
# ~ fig_timelapse.savefig(out_file, dpi=300)
# ~ fig_quiverzoom.savefig(out_file)

