import numpy as np
import gc
from scipy.interpolate import interpn
from scipy import ndimage

import pandas as pd
from scipy.interpolate import make_interp_spline


def div_fgradg(f,g): # generalize for f and g with dim 3
    
    
    out= np.zeros_like(g)
    
    filt = np.array([[0, 1, 0],[1, -4, 1], [0, 1, 0]])/2.   # careful that convolve reflects the filter w.r.t. what one would expect
    
    out += ndimage.convolve(f*g, filt, mode='nearest')
    
    
    filt = np.array([[0, 1, 0],[1, 0, 1], [0, 1, 0]])/2.    #
    
    out += f*ndimage.convolve(g, filt, mode='nearest')
    
    
    filt = np.array([[0, 1, 0],[1, 0, 1], [0, 1, 0]])/2.    # 
    
    out += -g*ndimage.convolve(f , filt, mode='nearest')
    
    return out
    
    

    

    
def derivatives_euler(pars, model, parsbact2=None): # implements the derivatives of a specified model
    
    #bacterial growth parameters
    r      =  pars['r']      # maximum growth rate of bacteria  ~ 1 = 40 mins doubling
    eps   =  pars['eps']     # biomass conversion factor
    k   =  pars['k']     # Monod constant of nutrients
    
    #phage growth parameters
    omega  = pars['omega']  # phage decay rate
    lambd  = pars['lambda']  #phage burst size
            
    #phage infection parameters
    phi= pars['phi']  # adsorption rate     
    
    Nstages =pars['Nstages'] # infection stages   
    eta= pars['eta']   # infection rate, inverse of latent time
    
    # motility parameters
    D_v=pars['D_v'] #diffusion coefficient for phage

    D_B=pars['D_B'] #Diffusion coefficient for bacteria
    am = pars['am'] # lowest sensed concentration
    ap = pars['ap'] # highest sensed concentration
    xi = pars['xi'] # bacteria chemosensing responsiveness
    D_R = pars["D_R"]   #diffusion coefficient for resources

    # attractant params, serine
    mu_1= pars["mu_1"]   # =0 
    mu_2= pars["mu_2"]      # max uptake rate
    k_attr= pars["k_attr"]  # monod constant 
    D_A= pars["D_A"]    #diffusion coefficient for atractants (equal to resources)
   
    # attractant params, aspartate
    mu_1a= pars["mu_1a"] # =0 
    mu_2a= pars["mu_2a"]   
    k_attra= pars["k_attra"]   
    D_A2= pars["D_A2"]  
    xi2  = pars["xi2"] # =xi 
    
    m_i = pars["mot_factor"]

    # chemotaxis parameters
    Delta_x = pars['Delta_x']
    
    # number of grid rows and cols
    Ncols = int(pars['L']/pars['Delta_x']) # number of grid rows and cols
    Nsites=int(Ncols**2)
    
    discr_thr= (10.**8.)/(Delta_x * Delta_x *0.5) # discreteness treshold for phage. The conversion factor transforms microns to mL
    
    def eta_ofR(R, eta_=eta, k_=k):
        # ~ return eta/3.
        return eta_*(R/(R+1*k_) ) #+ 0.3

    def beta_ofR(R, lambd_=lambd, k_=k):
        return 1 + lambd_*(R/(R+1*k_) )
    
    def discr_mask(same):
        return same**(10.)/(same**(10.) + discr_thr**(10.)  )

    def res_dyn(same, consumers, r_=r, eps_=eps, k_=k):
        return -r_ * eps_* same * consumers/(same + k_)  
        
    def attr_dyn(same, consumers, resources, k_attr_, mu_1_, mu_2_, r_=r, k_=k):
        return - (mu_1_*(r_ * resources/(resources + k_) ) + mu_2_)*  same * consumers/(same + k_attr_) 
        
    
    def bact_dyn(same, resources, phage_effect, r_=r, k_=k):
        return r_ * same * resources/(resources + k_)  - phage_effect*same
        
    def infection_stages_dyn(same, susceptibles, adv_rate, phage_effect):
        # ~ print "inf"
        dEi_dt= np.zeros_like(same)
        dEi_dt[0, :,:]= phage_effect*(susceptibles ) - adv_rate*same[0, :,:]
        dEi_dt[1:, :,:]= adv_rate*same[:-1, :,:] - adv_rate*same[1:, :,:]
        return dEi_dt
    
    def phage_dyn(same, phage_effect, burst_rate, infected, phage_effect2 =0, burst_rate_bact2=0, infected2=0):
        return burst_rate - phage_effect*(infected) - omega*same + burst_rate_bact2 - phage_effect2*(infected2)
   
    
    def diff(D, A):
        
        res=0
        stencil = np.array([[0.25, 0.5, 0.25],[0.5, -3, 0.5], [0.25, 0.5, 0.25]]) /(Delta_x**2.)  # careful that convolve reflects the filter w.r.t. what one would expect
        if (A.ndim==2):
            res= ndimage.convolve(D *A, stencil, mode='nearest')
        elif (A.ndim==3): # infected stages
            res = np.zeros_like(A)
            
            # there must be a more compact way using broadcasting, but it's not obvious how from the documentation of convolution packages
            for idx in range(A.shape[0]):
                res[idx,:,:] = ndimage.convolve(D *A[idx,:,:], stencil, mode='nearest')
         
        return res
        
        
    def chemotaxis_operator_Ping(B_or, R_or,xi_, am_=am, ap_=ap):
        
        B_ch = B_or # 
        R_ch = R_or
        
        res=0
        xiR = xi_
        if (B_ch.ndim==2):
            res= ( div_fgradg(xiR*B_ch,np.log((1+R_ch/am_)/(1+R_ch/ap_))) )/((Delta_x**2.))
        elif (B_ch.ndim==3): # infected stages
            res = np.zeros_like(B_ch)
            # there must be a more compact way using broadcasting, but it's not obvious how from the documentation of convolution packages
            for idx in range(B_ch.shape[0]):
                res[idx,:,:] = (  div_fgradg(xiR*B_ch[idx,:,:],np.log((1+R_ch/am_)/(1+R_ch/ap_))) )/((Delta_x**2.))
                
        else:
            print("wrong number of dimensions for bacteria array")
            sys.exit()
        
        return res
    
    def advance_one_step_attractant(dt, R,S,V,Ei,A, A2):
            
        etaR = eta_ofR(R)
        betaR = beta_ofR(R)
        discrhillV=discr_mask(V)
        
        B_tot = S +  np.sum(Ei, axis =0 ) # total number of bacteria
        
        
        
        Vc= pars["Pc"] # saturation threshold for phage binding 
        
        FV=V/(1+ V/Vc)
        if Vc <0:
            FV=V
            
        FS=1.
        FB=1.
        
                
        dRdt= res_dyn(R,B_tot) + diff(D_R,R)# + eps* Nstages*etaR*Ei[-1,:,:]*discrhillE[-1,:,:]/5.
        
        xieff=xi 
        dbact= D_B 
        dSdt= bact_dyn(S, R,  phi*discrhillV*FV*FS)   + diff(dbact,S) - chemotaxis_operator_Ping(S, A, xieff)- chemotaxis_operator_Ping(S, A2, xieff)
        D_v_eff=D_v

        dVdt= phage_dyn(discrhillV*V,  phi*discrhillV*FV, betaR*Nstages*etaR*Ei[-1,:,:], B_tot*FB) + diff(D_v_eff,V)

        dEidt= infection_stages_dyn(Ei, S, Nstages*etaR, phi*discrhillV*FV*FS)  + diff(dbact,Ei) - chemotaxis_operator_Ping(m_i*Ei, A, xi) - chemotaxis_operator_Ping(m_i*Ei, A2, xi2) 
        
        dAdt= attr_dyn(A,B_tot,R, k_attr, mu_1, mu_2) + diff(D_A,A) #+ mu_2* Nstages*etaR*Ei[-1,:,:]*discrhillE[-1,:,:]*  A /(A + k_attr) /2.
        dA2dt= attr_dyn(A2,B_tot,R, k_attra, mu_1a, mu_2a) + diff(D_A2,A2) #+ mu_2a* Nstages*etaR*Ei[-1,:,:]*discrhillE[-1,:,:]/5.

        R_new= R+ dt*dRdt 
        S_new= S+ dt*dSdt
        V_new= V+ dt*dVdt
        Ei_new= Ei+ dt*dEidt
        A_new= A+ dt*dAdt
        A2_new= A2+ dt*dA2dt

        
        return [R_new,S_new, V_new,Ei_new, A_new, A2_new]
        

        
    def advance_one_step_attractant_twobacts(dt, R,S,V,Ei,A, A2,S2, Ei2):
        
        r2          = parsbact2['r'] 
        eps2        = parsbact2['eps'] 
        k2          = parsbact2['k'] 
        lambd2  = parsbact2['lambda']  #burst size
                
        #phage infection parameters
        phi2= parsbact2['phi']  # adsorption rate     
        
        Nstages2 =Nstages # infection stages   
        eta2= parsbact2['eta']   # infection rate, inverse of latent time
        
    
        D_B2=parsbact2['D_B']
        am2 = parsbact2['am'] 
        ap2 = parsbact2['ap'] 
        xi_bact2 = parsbact2['xi'] 
    
        # attractant params
        mu_1_bact2= parsbact2["mu_1"]   
        mu_2_bact2= parsbact2["mu_2"]   
        k_attr_bact2= parsbact2["k_attr"]   
    
        mu_1a_bact2= parsbact2["mu_1a"]   
        mu_2a_bact2= parsbact2["mu_2a"]   
        k_attra_bact2= parsbact2["k_attra"]   
        
        m_i2 = parsbact2["mot_factor"]
   
        etaR  = eta_ofR(R)
        betaR = beta_ofR(R)
        etaR2 = eta_ofR(R, eta2, k2)
        betaR2 = beta_ofR(R, lambd2, k2)
        discrhillV=discr_mask(V)
        
        
        
        B_tot = S +  np.sum(Ei, axis =0 ) +  np.sum(Ei2, axis =0 ) + S2 # total number of bacteria

        
        Vc= pars["Pc"] # saturation threshold for phage binding 
        
        FV=V/(1+ V/Vc)
        if Vc <0:
            FV=V
        FS=1.
        FB=1.
        
                
        dRdt= res_dyn(R,S +  np.sum(Ei, axis =0 )) +res_dyn(R, np.sum(Ei2, axis =0 ) + S2, r2, eps2, k2) + diff(D_R,R)# + eps* Nstages*etaR*Ei[-1,:,:]*discrhillE[-1,:,:]/5.
        
        xieff=xi 
        dbact= D_B 
        dSdt= bact_dyn(S, R,  phi*discrhillV*FV*FS)   + diff(dbact,S) - chemotaxis_operator_Ping(S, A, xieff)- chemotaxis_operator_Ping(S, A2, xieff)
        dS2dt= bact_dyn(S2, R, phi2*discrhillV*FV, r2, k2)   + diff(D_B2,S2) - chemotaxis_operator_Ping(S2, A, xi_bact2, am2, ap2)- chemotaxis_operator_Ping(S2, A2, xi_bact2, am2, ap2)
        
       
        D_v_eff=D_v

        dVdt= phage_dyn(discrhillV*V,  phi*discrhillV*FV, betaR*Nstages*etaR*Ei[-1,:,:], (S +  np.sum(Ei, axis =0 )), phi2*discrhillV*FV, betaR2*Nstages2*etaR2*Ei2[-1,:,:], (S2 +  np.sum(Ei2, axis =0 ))) + diff(D_v_eff,V)


        # ~ print (np.max(Ei),np.max(Ei[-1,:,:]), np.max(etaR),np.max(phi*discrhillV*FV*FS), phi)
        dEidt= infection_stages_dyn(Ei, S, Nstages*etaR, phi*discrhillV*FV*FS)  + diff(dbact,Ei) - chemotaxis_operator_Ping(m_i*Ei, A, xi) - chemotaxis_operator_Ping(m_i*Ei, A2, xi) 
        
        dEi2dt= infection_stages_dyn(Ei2, S2, Nstages2*etaR2, phi2*discrhillV*FV*FS)  + diff(D_B2,Ei2) - chemotaxis_operator_Ping(m_i2*Ei2, A, xi_bact2, am2, ap2) - chemotaxis_operator_Ping(m_i2*Ei2, A2, xi_bact2, am2, ap2) 
        
        dAdt= attr_dyn(A,S +  np.sum(Ei, axis =0 ), R, k_attr, mu_1, mu_2) + attr_dyn(A,np.sum(Ei2, axis =0 ) + S2,R, k_attr_bact2, mu_1_bact2, mu_2_bact2, r2, k2) + diff(D_A,A) 
        dA2dt= attr_dyn(A2,S +  np.sum(Ei, axis =0),R, k_attra, mu_1a, mu_2a)+ attr_dyn(A2,np.sum(Ei2, axis =0 ) + S2,R, k_attra_bact2, mu_1a_bact2, mu_2a_bact2, r2, k2) + diff(D_A2,A2) 
        
        R_new= R+ dt*dRdt 
        S_new= S+ dt*dSdt
        V_new= V+ dt*dVdt
        Ei_new= Ei+ dt*dEidt
        A_new= A+ dt*dAdt
        A2_new= A2+ dt*dA2dt
        S2_new= S2+ dt*dS2dt
        Ei2_new= Ei2+ dt*dEi2dt

        
        
        return [R_new,S_new, V_new,Ei_new, A_new, A2_new, S2_new, Ei2_new]
        
        
    def return_derivatives_and_fluxes( R,S,V,E,A, A2): # assuming Ei is uniform
        
        etaR = eta_ofR(R)
        betaR = beta_ofR(R)
        discrhillV=discr_mask(V)
        
        
        B_tot = S +  E   # total number of bacteria
        Vc= pars["Pc"] # saturation threshold for phage binding 
        
        FV=V/(1+ V/Vc)
        if Vc <0:
            FV=V
        FS=1.
        FB=1.
        
                
        dRdt= res_dyn(R,B_tot) + diff(D_R,R)# + eps* Nstages*etaR*Ei[-1,:,:]*discrhillE[-1,:,:]/5.
        
        xieff=xi 
        dbact= D_B 
        dSdt= bact_dyn(S, R,  phi*discrhillV*FV*FS)   + diff(dbact,S) - chemotaxis_operator_Ping(S, A, xieff)- chemotaxis_operator_Ping(S, A2, xieff)
        
        D_v_eff=D_v

        dVdt= phage_dyn(discrhillV*V,  phi*discrhillV*FV, betaR*etaR*E, B_tot*FB) + diff(D_v_eff,V)

        dAdt= attr_dyn(A,B_tot,R, k_attr, mu_1, mu_2) + diff(D_A,A) #+ mu_2* Nstages*etaR*Ei[-1,:,:]*discrhillE[-1,:,:]*  A /(A + k_attr) /2.
        dA2dt= attr_dyn(A2,B_tot,R, k_attra, mu_1a, mu_2a) + diff(D_A2,A2) #+ mu_2a* Nstages*etaR*Ei[-1,:,:]*discrhillE[-1,:,:]/5.
        
        return [dRdt,dSdt, dVdt, dAdt, dA2dt, chemotaxis_operator_Ping(m_i*B_tot, A2, xi2), diff(dbact,B_tot), r * S * R/(R + k), phi*discrhillV*FV*FS*S, etaR*E,  np.log((1+A2/am)/(1+A2/ap))] # returns derivatives, plus chemotaxis 2, diffusion rate, bacteria growth rate, phage infection rate, cell lysis rate, attractant 2 sensing field
        
        
    if model == 'attractant':
        return advance_one_step_attractant
    elif model == 'twobacts':
        return advance_one_step_attractant_twobacts
    elif model == 'derivatives':
        return return_derivatives_and_fluxes
    else:
        print("model not implemented")
        sys.exit()


def bact_ext_prep(pars):
    

    Nstages =pars['Nstages'] # infection stages   
    Ncols = int(pars['L']/pars['Delta_x']) # number of grid rows and cols
    Nsites=int(Ncols**2)    
    
    # extinction on the whole grid, then stop the integration. I need to implement a "local extinction" rule
    def bact_ext(t, z): return np.sum(z[Nsites:2*Nsites]) + np.sum(z[3*Nsites:(Nstages+3)*Nsites])  - 1e4 # threshold at 1e4. if Delta_x 100 um, area is 1e4 um**2, vs 1e8 um**2 units (cm**2)
    
    return bact_ext



    
def bact_ext_inc(pars):
    Nstages =pars['Nstages'] # infection stages   
    def bact_ext(t, z): 
        R, S, V  = z[:3]
        Ei= z[3:3+Nstages]
        B_tot = S +  np.sum(Ei, axis =0 ) # total number of bacteria
        return B_tot - 1. # threshold at 10 /mL. If the lung volume is 100 ul that corresponds to 1 cell. Perhaps it's 400 uL

    return bact_ext
    
    
    
def bact_ext_lambda(t, z, pars):
    Nstages =pars['Nstages'] # infection stages           
            
    R, S, V  = z[:3]
    Ei= z[3:3+Nstages]
    B_tot = S +  np.sum(Ei, axis =0 ) # total number of bacteria
    
    return B_tot - 1. # threshold at 10 /mL. If the lung volume is 100 ul that corresponds to 1 cell. Perhaps it's 400 uL


def coarsen2D(original_array, c):
    temp = original_array.reshape((int(original_array.shape[0] // c), c, int(original_array.shape[1] // c), c))
    coarse_arr = np.mean(temp, axis=(1,3))
    return coarse_arr
  

def gridint(x, y, z, xi, yi):
    "Convert 3 column data to matplotlib grid"
    
    print((min(xi), min(yi)))
    print((max(xi), max(yi), max(y)))
    Z = interpn((x, y), z, (xi[None,:], yi[:,None]), method='linear')
    return Z
    
    
    
def rotate(origin, px, py, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    ox, oy = origin
    # ~ px, py = point

    qx = ox + np.cos(angle) * (px - ox) - np.sin(angle) * (py - oy)
    qy = oy + np.sin(angle) * (px - ox) + np.cos(angle) * (py - oy)
    return qx, qy
    
