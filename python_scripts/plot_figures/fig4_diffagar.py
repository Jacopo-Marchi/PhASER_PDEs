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


fontprops = fm.FontProperties(size=28)

dir_fig='../../sims_for_fig4_new' # directory with input files


df = pd.read_excel('{inp}/../growth_data/Scanner_OD.xlsx'.format(inp=dir_fig), skiprows=1, usecols='A,B',  header = None,  engine='openpyxl', sheet_name="CFU Calibration")


pxtoCFUs=df.to_numpy(copy=True)[:8,:]

pxtoCFUs=pxtoCFUs[np.argsort(pxtoCFUs[:,1]),:]

print(pxtoCFUs)
print(pxtoCFUs[0,:])

# ~ f = interpolate.interp1d(pxtoCFUs[:,1], pxtoCFUs[:,0])
f = make_interp_spline(pxtoCFUs[:,1], pxtoCFUs[:,0], k=3)

############ PLOT FRAMES

dirstrings=['{inp}/F8_t4_Vc1e7_agar0d2','{inp}/F8_t7_Vcinf_agar0d2','{inp}/F8_t4_Vc1e7_agar0d3','{inp}/F8_t7_Vcinf_agar0d3']
# ~ dirstrings=['{inp}/t4_F8_xi152e5_DB_32e4','{inp}/t7_F8_xi152e5_DB_32e4','{inp}/t4_F8_xi147e5_DB_32e4','{inp}/t7_F8_xi147e5_DB_32e4']
t_ints=[38,34,64,58]
# ~ t_ints=[20,20,20,20]

    
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
    
    t_estspeed=[22]
    t_estspeedm2=[tm2-2 for tm2 in t_estspeed]
    t_int=t_ints[idir]
    speedest=0
    wavespeeds=np.zeros_like(times)


    
    
    for it, t in enumerate(times):
        
        # ~ if int(t) == t_int or it == 0 or int(t) in t_estspeed + t_estspeedm2:

        print (t_int, int(t))
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

            
            print(wavepos)
            print(frontpos)
            
            if waveposprec is None:
                waveposprec=wavepos
            if frontposprec is None:
                frontposprec=frontpos
                
            # ~ wavespeed= (wavepos - waveposprec)/(times[it]-times[it-1])
            wavespeed= (frontpos - frontposprec)/(times[it]-times[it-1])
            print("time: ",time)
            print("wavespeed: ",wavespeed)
              
            if wavespeed>0 and int(t) in t_estspeed:
                speedest=wavespeed
            
            
            if wavespeed>0:
                wavespeeds[it]=wavespeed
        
        
        bactprec=B_tot
        intprec=intensities
        waveposprec=wavepos
        frontposprec=frontpos
    

        if int(t) == t_int:
            
            print ("PLOTTING ")
            print (t_int, int(t))
            
            fig = plt.figure(figsize=(4,4))
            grid = gridspec.GridSpec(1, 1, left=0.03, right=0.97, top=0.95, bottom=0.03,
                 wspace=0.4, hspace=0.35)
            labeled_axes = []
            ax = plt.Subplot(fig, grid[0, 0])
            fig.add_subplot(ax)
            
            
            if '0d3' in dirstr and L > 120000:
                ax.set_xlim(xmin=-100, xmax= 20)
                ax.set_ylim(ymin=-100, ymax= 20)
                # ~ ax.annotate("{:.1f} mm/hr".format(speedest), xy=(-100+10*(12/20), 20- 15*(12/20)), xycoords='data',xytext=(0, 5), textcoords='offset pixels', color='k', fontsize=16)   
                
            # ~ elif np.amin(xx) < -60:
                # ~ ax.annotate("{:.1f} mm/hr".format(speedest), xy=(-90, 85), xycoords='data',xytext=(0, 5), textcoords='offset pixels', color='k', fontsize=16)  
                          
            # ~ else:
                # ~ ax.annotate("{:.1f} mm/hr".format(speedest), xy=(-60+10*(12/20), 60- 15*(12/20)), xycoords='data',xytext=(0, 5), textcoords='offset pixels', color='k', fontsize=16)             
    

                
            v= plot_timelapse_intensity(xx,yy,intensities,ax,0) #,rescdim=10/6
            
            ax.plot((21000  -L/2.)/1000, (21000  -L/2.)/1000, marker='x', markersize=12, color='k', linestyle='')
            ax.plot((6000  -L/2.)/1000, (6000  -L/2.)/1000, marker='o',  markerfacecolor='none', markersize=12,   color='k', linestyle='')
            # ~ ax.annotate("Front speed {:.1f} mm/hr".format(speedest), xy=(-90, 85), xycoords='data',xytext=(0, 5), textcoords='offset pixels', color='k', fontsize=16)             
            
            out_file='{out}/last_intensity.png'.format(out=dir_out_plots_tot)
            #    print out_file
            fig.savefig(out_file, dpi=300)
            fig.clf()
            plt.close('all')
        
        
        
                
    
    
    fig = plt.figure(figsize=(4,4))
    grid = gridspec.GridSpec(1, 1, left=0.15, right=0.7, top=0.91, bottom=0.22, wspace=0.4, hspace=0.35)
    labeled_axes = []
        
    ax = plt.Subplot(fig, grid[0, 0])
    fig.add_subplot(ax)
    labeled_axes.append(ax)
    
    ax.plot(times, wavespeeds , linestyle='-', color='r', label='R')
    
    ax.set_xlabel('time (h)')
    ax.set_ylabel('wave speed')
    ax.xaxis.labelpad = axis_labelpad
    ax.yaxis.labelpad = axis_labelpad
    ax.legend(frameon=False, ncol=1, columnspacing=0.5, handletextpad=0.2,
        loc='upper right', bbox_to_anchor=(1.5, 1.18))
    
    mpsetup.despine(ax) 
        
    out_file='{out}/wavespeed.pdf'.format(out=dir_out_plots_tot)
    fig.savefig(out_file) 
    
 
        
    data_fin = np.array([times, wavespeeds ])
     
    file_out='{inp}/wavespeeds.txt'.format(inp=dir_out_plots_tot) # saves speeds, the reference value is taken at 22 hours
    
    with open(file_out,'w') as f_handle:
        np.savetxt(f_handle, data_fin.T, fmt='%40.15f',  header="x \t wave speed")
        f_handle.write("\n")
        
     
    
