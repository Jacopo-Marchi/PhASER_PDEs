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



# ============================================================
# MODEL FUNCTIONS (Choua & Bonachela 2019, eqs. 1-2)
# ============================================================
def eclipse(mu_norm, E_inf, E0, alpha_E):
    """E(mu) = E_inf + E0 * exp(-alpha_E * mu/mu_max)  [minutes]"""
    return E_inf + E0 * np.exp(-alpha_E * mu_norm)

def maturation(mu_norm, M_inf, M0, alpha_M):
    """M(mu) = M_inf / (1 + exp(-alpha_M*(mu/mu_max - M0)))  [virions/min]"""
    return M_inf / (1.0 + np.exp(-alpha_M * (mu_norm - M0)))


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



dirstrings=['{inp}/t7_F9_lowK','{inp}/t4_F9_lowK','{inp}/t7_F9_lowK_plasticvirus','{inp}/t4_F9_lowK_plasticvirus']
t_ints=[26,32,26,32]
titstrings=[r'Phage T7, linear',r'Phage T4, linear',r'Phage T7, nonlinear',r'Phage T4, nonlinear']


figallee, axs_allee = plt.subplots(nrows=3, ncols=2, figsize=(8,12),  layout='constrained')# , layout='constrained',


    
labeled_axes = []
labeled_axes.append(axs_allee[0,0])
# ~ labeled_axes.append(axs_fig2[1,0])
labeled_axes.append(axs_allee[0,1])
labeled_axes.append(axs_allee[1,0])
labeled_axes.append(axs_allee[1,1])
labeled_axes.append(axs_allee[2,0])
labeled_axes.append(axs_allee[2,1])


 
mu_norm = np.linspace(0.01, 1.0, 300)


# ~ MU_MAX_MODEL = 1.6       # h^-1, reference mu_max for model lines
# ~ FRAC_SHORT   = 3.5 / (3.5 + 5.0)   # = 1/3, x-range factor for "short" variant

dir_io=dirstrings[0].format(inp=dir_fig)  # directory with input files

file_params='{inp}/parameters.json'.format(inp=dir_io)
pars_T7 = {}
with open(file_params,'r') as f_handle:
    pars_T7 = json.loads(f_handle.read())
    
dir_io=dirstrings[1].format(inp=dir_fig)  # directory with input files

file_params='{inp}/parameters.json'.format(inp=dir_io)
pars_T4 = {}
with open(file_params,'r') as f_handle:
    pars_T4 = json.loads(f_handle.read())
    
k_ =     pars_T7['k'] # directory with input files
MU_MAX_MODEL =     pars_T7['r'] # directory with input files
print(k_, MU_MAX_MODEL)


datasets = {
    'Phage T7, nonlinear': {
        'burst': pars_T7['lambda'],      # 
        'eta': pars_T7['eta'],      # 
        # Eclipse
        'E_inf': 15.33, 'E0': 82.71, 'alpha_E': 5.9, # Inferred from You
        'M_inf': 6.7, 'M0': 0.64, 'alpha_M': 9.70,
        'color': 'r',  'ls': '-'
    },
    'Phage T4, nonlinear': { 
        'burst': pars_T4['lambda'],      # 
        'eta': pars_T4['eta'],      # 
        'E_inf': 18.26, 'E0': 94.18, 'alpha_E': 7.76, # Inferred from Hadas and Golec datasets
        'M_inf': 47.32, 'M0': 0.65, 'alpha_M': 12.03,
        'color': 'b', 'ls': '-'
}
}       

FRAC_SHORT   = 3500 / (3500 + k_)  
print(FRAC_SHORT)

# ~ sys.exit()

# ---- Choua & Bonachela curves ----
for label, d in datasets.items():
    col    = d['color']
    ls     = d['ls']
    lambd_ = d['burst']
    eta_ = d['eta']
    curve_label = label

        
    E_star = eclipse(FRAC_SHORT, d['E_inf'], d['E0'], d['alpha_E'])
    L_model = 60.0 / (eta_* FRAC_SHORT)          # min
    w_norm = 1.0 / (L_model - E_star) # min^-1
    E_curve   = eclipse(mu_norm, d['E_inf'], d['E0'], d['alpha_E'])
    L_curve   = 1.0 / w_norm + E_curve   # [min]
    eta_curve = 60.0 / L_curve    # [h^-1]
    
    M_star = maturation(FRAC_SHORT,d['M_inf'], d['M0'], d['alpha_M'])
    B_model = 1 + lambd_* FRAC_SHORT
    w_norm   = M_star / B_model   
    M_curve   = maturation(mu_norm, d['M_inf'], d['M0'], d['alpha_M'])
    B_curve   = M_curve / w_norm
    
    
    axs_allee[0,0].plot(mu_norm, eta_curve, ls=ls, color=col, lw=2, label=curve_label)
    axs_allee[0,1].plot(mu_norm,   B_curve,   ls=ls, color=col, lw=2, label=curve_label)
 
 
 
    mu_lin = mu_norm * MU_MAX_MODEL   # actual mu in h^-1
 
    eta_lin = (eta_  / MU_MAX_MODEL) * mu_lin
    B_lin   = 1.0 + (lambd_/ MU_MAX_MODEL) * mu_lin
 
    axs_allee[0,0].plot(mu_norm, eta_lin, ls='--', color=col,
                lw=2, label=curve_label.replace('nonlinear', 'linear'))
 
    axs_allee[0,1].plot(mu_norm, B_lin, ls='--', color=col,
              lw=2, label=curve_label.replace('nonlinear', 'linear'))
 
    # ---- formatting ----
    xlabel = r'Normalized growth rate'
 
    axs_allee[0,0].set_xlabel(xlabel, fontsize=16)
    axs_allee[0,0].set_ylabel(r'Lysis rate (h$^{-1}$)', fontsize=16)
    axs_allee[0,0].legend(frameon=False, fontsize=12, loc='upper left')
    axs_allee[0,0].set_xlim(0, FRAC_SHORT)
    axs_allee[0,0].set_ylim(0, 4.9)
    # ~ axs_allee[0,0].set_xticklabels(['0' if t == 0.0 else f'{t:.1f}' for t in axs_allee[0,0].get_xticks()])
    axs_allee[0,0].xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, _: '0' if x == 0.0 else f'{x:.1f}'))
 
    axs_allee[0,1].set_xlabel(xlabel, fontsize=16)
    axs_allee[0,1].set_ylabel('Burst size (PFU/CFU)', fontsize=16)
    axs_allee[0,1].legend(frameon=False, fontsize=12, loc='upper left')
    axs_allee[0,1].set_xlim(0, FRAC_SHORT)
    axs_allee[0,1].set_ylim(0, 69)
    # ~ axs_allee[0,1].set_xticklabels(['0' if t == 0.0 else f'{t:.1f}' for t in axs_allee[0,1].get_xticks()])
    axs_allee[0,1].xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, _: '0' if x == 0.0 else f'{x:.1f}'))


print (mu_norm)

# ~ sys.exit()
    
for idir, dirstr in enumerate(dirstrings):

    dir_io=dirstr.format(inp=dir_fig)  # directory with input files
    
    dir_in_tot='{inp}/data'.format(inp=dir_io) # directory with output plots
    dir_in_frames='{inp}/frames'.format(inp=dir_in_tot) # directory with output plots
    dir_out_plots_tot='{inp}/plots'.format(inp=dir_io) # directory with output plots
        
        
    
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
        
    Ncols  = int(pars['L']/pars['Delta_x']) # number of grid rows and cols
    Nrows  = int(pars['Ly']/pars['Delta_x']) # number of grid rows and cols
    Nsites = int(Ncols*Nrows)
    
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
    t_int=t_ints[idir]
    
    speedest=0
    wavespeeds=np.zeros_like(times)
    
    
    for it, t in enumerate(times):
        
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
            
            
            
            xlin=xx[id_first_bacts,id_first_bacts_x:]

            
            
            wavepos=xlin[np.argmax(wave[id_first_bacts,id_first_bacts_x:])]
            frontpos=xlin[np.argmax(B_tot[id_first_bacts,id_first_bacts_x:]<1e7)]

            
            if waveposprec is None:
                waveposprec=wavepos
            if frontposprec is None:
                frontposprec=frontpos
                
            # ~ wavespeed= (wavepos - waveposprec)/(times[it]-times[it-1])
            wavespeed= (frontpos - frontposprec)/(times[it]-times[it-1])
            print("time: ",time)
            print("wavespeed: ",wavespeed)
              
            if wavespeed>0:
                wavespeeds[it]=wavespeed
        
        
        bactprec=B_tot
        intprec=intensities
        waveposprec=wavepos
        frontposprec=frontpos
    
    
        if int(t) == t_int:
            
            row= int(idir/2) + 1
            col= np.mod(idir,2)
                
            # ~ heatmap = ax.pcolormesh(xx, yy, intensities, cmap='gray',vmin=maxbacts_convert, vmax=no_bactcontrol) # faster than pcolor
            heatmap = axs_allee[row,col].pcolormesh(xx, yy, intensities, cmap='gray',vmin=0, vmax=255, shading='auto', linewidth=0,rasterized=True) # faster than pcolor
            # ~ ax.set_xlabel(r'x (mm)')
            # ~ ax.set_ylabel(r'y (mm)')
            
            axs_allee[row,col].plot((21000  -L/2.)/1000, (21000  -L/2.)/1000, marker='x', markersize=6, color='k', linestyle='')
            axs_allee[row,col].plot((6000  -L/2.)/1000, (6000  -L/2.)/1000, marker='o',  markerfacecolor='none', markersize=6,   color='k', linestyle='')
            
          
            
            # ~ cbar1 = fig.colorbar(heatmap,ax=ax)
            
            # ~ cbar1.ax.get_yaxis().labelpad = 15
            # ~ cbar1.ax.set_ylabel('pixel intensity', rotation=270)
                
              
            heatmap = axs_allee[row,col].pcolormesh(xx, yy,   wavesum , cmap='gray_r', alpha=0.5) # faster than pcolor
            
            axs_allee[row,col].set_aspect('equal')
    
            scalebar = AnchoredSizeBar(axs_allee[row,col].transData,
                                   25, '25 mm', 'lower right', 
                                   pad=0.2,
                                   color='black',
                                   frameon=False,
                                   size_vertical=0.3, fontproperties=fontprops)
            axs_allee[row,col].add_artist(scalebar)
          
            trans = transforms.blended_transform_factory(axs_allee[row,col].transData, axs_allee[row,col].get_xticklabels()[0].get_transform())
            
            axs_allee[row,col].annotate(titstrings[idir], xy=(-90, 85), xycoords='data',xytext=(0, 5), textcoords='offset pixels', color='k', fontsize=16)

    
            # ~ ax.annotate(str(int(t))+" hr", xy=(65, 88), xycoords='data',xytext=(0, 5), textcoords='offset pixels', color='k', fontsize=12)    #
        #
            axs_allee[row,col].tick_params(left = False, right = False , labelleft = False ,
                        labelbottom = False, bottom = False)
                        
                
        
            
        data_fin = np.array([times, wavespeeds ])
         
        file_out='{inp}/wavespeeds.txt'.format(inp=dir_out_plots_tot) # saves speeds, the reference value is taken at 22 hours
        
        with open(file_out,'w') as f_handle:
            np.savetxt(f_handle, data_fin.T, fmt='%40.15f',  header="x \t wave speed")
            f_handle.write("\n")
                      
#### finish figure ####
labeldict = dict(labelstyle=r'%s', fontsize=26,
     xycoords=('axes fraction'), fontweight = 'bold')
mpsetup.label_axes([labeled_axes[0]], labels='A', xy=(-0.1, 0.95), **labeldict)
mpsetup.label_axes([labeled_axes[1]], labels='B', xy=(-0.1, 0.95), **labeldict)
mpsetup.label_axes([labeled_axes[2]], labels='C', xy=(-0.1, 0.95), **labeldict)
mpsetup.label_axes([labeled_axes[3]], labels='D', xy=(-0.1, 0.95), **labeldict)
mpsetup.label_axes([labeled_axes[4]], labels='E', xy=(-0.1, 0.95), **labeldict)
mpsetup.label_axes([labeled_axes[5]], labels='F', xy=(-0.1, 0.95), **labeldict)
out_file='{out}/SIfig_plasticvir_lowK.png'.format(out=dir_fig)
#    print out_file
figallee.savefig(out_file, dpi=300)
figallee.clf()
plt.close('all')
    
 
    
 
