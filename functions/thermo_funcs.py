###########################################################################
#                                                                         #
#    Courtesy of Jake Seeley, Nov 2016.                                   #
#                                                                         #
###########################################################################


import numpy as np
import scipy.optimize as opt 
from thermo_constants import *

##  Saturation vapor pressure with respect to liquid  ###################
def pvstar_l(T):
    return (e_0*(T/T_0)**((c_pv - c_vl)/R_v))* \
           np.exp(((E_0v - T_0*(c_vv - c_vl))/R_v)* \
           (1./T_0 - 1./T))

##  Saturation vapor pressure with respect to solid  ####################
def pvstar_s(T,icemode=1):
    if icemode==1:
        c_vs2 = c_vs
        E_0s2 = E_0s
        s_0s2 = s_0s
    elif icemode==0:
        c_vs2 = c_vl
        E_0s2 = 0.
        s_0s2 = 0. 
    return (e_0*(T/T_0)**((c_pv - c_vs2)/R_v))* \
           np.exp(((E_0v + E_0s2 - T_0*(c_vv - c_vs2))/R_v)* \
           (1./T_0 - 1./T))

##  Saturation vapor pressure with transition at triple point temp  #####
def pvstar(T,icemode=1):
    if np.isscalar(T):
        if T > T_0:
            return pvstar_l(T)
        else:
            return pvstar_s(T,icemode=icemode)
    elif isinstance(T,np.ndarray):
        shape = T.shape
        pvstar_array = np.nan*np.empty(shape)
        pvstar_array[T>T_0]  = pvstar_l(T[T>T_0])
        pvstar_array[T<=T_0] = pvstar_s(T[T<=T_0],icemode=icemode) 
        return pvstar_array
    else:
        print 'Error: supplied temperature data not array or scalar'
        return
        
##  Saturation specific humidity  #######################################
def qvstar(p,T,qt=0.,ql=0.,qs=0.,liquid=False,solid=False,icemode=1):
    if (np.any(ql < 0.) or np.any(qs < 0) or np.any(qt < 0.)):
        print 'Error: negative mass fractions supplied'
        print 'qt: ',qt
        print 'ql: ',ql
        print 'qs: ',qs
        return
    if icemode==1:
        c_vs2 = c_vs
        E_0s2 = E_0s
        s_0s2 = s_0s
    elif icemode==0:
        c_vs2 = c_vl
        E_0s2 = 0.
        s_0s2 = 0. 
    if liquid:
        dcv = c_vv - c_vl
        dcp = c_pv - c_vl
        E0  = E_0v
    elif solid:
        dcv = c_vv - c_vs2
        dcp = c_pv - c_vs2
        E0  = E_0v + E_0s2
    else:
        warm = 0.5 + 0.5*np.sign(T-T_0)
        dcv  = (c_vv - c_vs2) + warm*((c_vv - c_vl) - (c_vv - c_vs2))
        dcp  = (c_pv - c_vs2) + warm*((c_pv - c_vl) - (c_pv - c_vs2))
        E0   = (E_0v + E_0s2) - warm*E_0s2
    if (qt != 0.):
        return (1. - qt)*(R_d/R_v)/(p/e_0 * (T_0/T)**(dcp/R_v) * \
             np.exp(-(E0 - dcv*T_0)/R_v*(1./T_0 - 1./T)) - 1.)
    else:
        ql2 = ql
        qs2 = qs
        return (1. - ql2 - qs2)/(R_v*p/(R_d*e_0)*(T_0/T)**(dcp/R_v)* \
               np.exp(-(E0 - dcv*T_0)/R_v*(1./T_0 - 1./T)) - eps)
               
def qvstar_mixedphase(p,T,qt,icemode=1):
    eta = mixedphase_eta(T)
    return (1.-eta)*qvstar(p,T,qt=qt,liquid=True,icemode=icemode)+eta*qvstar(p,T,qt=qt,solid=True,icemode=icemode)

def mixedphase_eta(T,Tfreeze=233.15):
    if np.isscalar(T):
        if (T <= Tfreeze):
            eta = 1.
        elif (Tfreeze < T and T < T_0):
            eta = (T_0-T)/(T_0 - Tfreeze)
        else:
            eta = 0.
        return eta
    elif isinstance(T,np.ndarray):
        shape = len(T)
        eta_array = np.nan*np.empty(shape)
        for i in range(shape):
            if T[i] <= Tfreeze:
                eta_array[i] = 1.
            elif (Tfreeze < T[i] and T[i] < T_0):
                eta_array[i]  = (T_0-T[i])/(T_0-Tfreeze)
            else:
                eta_array[i] = 0.     
        return eta_array
    else:
        print 'Error: supplied temperature data not array or scalar'
        return
    
def MSE(T,z,qv,ql,qs,icemode=1):
    if icemode==1:
        c_vs2 = c_vs
        E_0s2 = E_0s
    elif icemode==0:
        c_vs2 = c_vl
        E_0s2 = 0. 
    qa = 1. - qv - ql - qs
    cpm = qa*c_pd + qv*c_pv + ql*c_vl + qs*c_vs2
    return cpm*(T-T_0) + qv*(E_0v + R_v*T_0) - qs*E_0s2 + gg*z
    
def MSE_simple(z,T,qv):
    return gg*z + c_pd*T + L_v*qv
        
##  Equivalent potential temperature  ###################################
def theta_e(p,T,qv,ql=0.,qs=0.,p_0=p_0,icemode=1):
    qa = 1. - qv - ql - qs
    rv = qv/qa
    rl = ql/qa
    rs = qs/qa
    pa = (R_d/(R_d + rv*R_v))*p
    pv = rv*(R_v/(R_d + rv*R_v))*p
    if (icemode == 1):
        c_vs2 = c_vs
        E_0s2 = E_0s
        s_0s2 = s_0s
    elif (icemode == 0):
        c_vs2 = c_vl
        E_0s2 = 0.
        s_0s2 = 0.
    if np.isscalar(pv):
    	if (pv == 0.):
    		vap_pres_ratio = 0.
    	else:
        	vap_pres_ratio = e_0/pv
    elif isinstance(pv,np.ndarray):
    	shape = pv.shape
    	vap_pres_ratio = np.nan*np.empty(shape)
    	vap_pres_ratio[pv==0.] = 0.
    	vap_pres_ratio[pv!=0.] = e_0/pv
    	
    theta_e_vals = T*(p_0/pa)**(R_d/c_pd)*(T/T_0)** \
                  ((rv*c_pv + rl*c_vl + rs*c_vs2)/c_pd)* \
                  (vap_pres_ratio)**(rv*R_v/c_pd)*np.exp((rv*s_0v - rs*s_0s2)/c_pd)
    return theta_e_vals
    
def fallout(fallout_factor,p,T,qv_series,ql_series,qs_series,icemode=1,verbose=False):
    thetae_init = theta_e(p,T,qv_series[1],ql_series[1],qs_series[1],icemode=icemode)
    qt_init = qv_series[1]+ql_series[1]+qs_series[1]
    qt_min = qvstar_mixedphase(p,T,qt=qt_init,icemode=icemode)
    dqv = qv_series[1]-qv_series[0]
    dql = ql_series[1]-ql_series[0]
    dqs = qs_series[1]-qs_series[0]
    dqcon = dql + dqs
    
    if verbose:
        print 'theta_e before fallout: '+str(thetae_init)
        print 'qt before fallout: '+str(qt_init)
        print 'Minimum qt: '+str(qt_min)
        print 'dqv: '+str(dqv)
        print 'dql: '+str(dql)
        print 'dqs: '+str(dqs)
        print 'new condensates: '+str(dqcon)
    # If no new condensates, no fallout
    if (dqcon <= 0.):
        if verbose:
            print 'No new condensates.'
        ql_new = ql_series[1]
        qs_new = qs_series[1]
        thetae_new = thetae_init
        qt_new = qt_init
    else:
        # more liquid and ice
        if (dql >= 0. and dqs >= 0.):
            dql_new = (1. - fallout_factor)*dql
            dqs_new = (1. - fallout_factor)*dqs
        # less liquid, more ice
        elif (dql < 0. and dqs > 0.):
            dql_new = dql
            dqs_new = dqs - fallout_factor*dqcon
        # more liquid, less ice
        elif (dql > 0. and dqs < 0.):   
            dql_new = dql - fallout_factor*dqcon
            dqs_new = dqs_new
        else:
            print 'warning: other case'
        qv_new = qv_series[1]
        ql_new = ql_series[0] + dql_new
        qs_new = qs_series[0] + dqs_new
        qt_new = qv_series[1] + ql_new + qs_new
        dqcon_new = ql_new - ql_series[0] + qs_new - qs_series[0]
        if (qt_new < qt_min):
            print 'warning: parcel now subsaturated'
            qt_new = qt_min
            qv_new = qt_min
            ql_new = 0.
            qs_new = 0.
        thetae_new = theta_e(p=p,T=T,qv=qv_new,ql=ql_new,qs=qs_new)
        if verbose:
            print 'New ql: '+str(ql_new)
            print 'New qs: '+str(qs_new)
            print 'New thetae: '+str(thetae_new)
            print 'New qt: '+str(qt_new)
            print 'Fraction fallout: '+str(1. - dqcon_new/dqcon)
            
    return_dict = {}
    return_dict['ql_new']     = ql_new
    return_dict['qs_new']     = qs_new
    return_dict['thetae_new'] = thetae_new
    return_dict['qt_new']     = qt_new
        
    return (return_dict) 
    
def partition_mixedphase(p,T,qt,icemode=1):
    # print "partition_mixedphase: qt =",qt,"qvstar_mixedphase =",qvstar_mixedphase(p,T,qt=qt,icemode=icemode)
    if (qt >= qvstar_mixedphase(p,T,qt=qt,icemode=icemode)):
    	if icemode==1:
        	qv = qvstar_mixedphase(p,T,qt=qt,icemode=icemode)
        	eta = mixedphase_eta(T)
        	cond = qt-qv
        	ql = (1.-eta)*cond
        	qs = eta*cond
        elif icemode==0:
        	qv = qvstar_mixedphase(p,T,qt=qt,icemode=icemode)
        	ql = qt - qv
        	qs = 0.
    else:
        qv = qt
        ql = 0.
        qs = 0.
    return_dict = {}
    return_dict['qv'] = qv
    return_dict['ql'] = ql
    return_dict['qs'] = qs
    return return_dict
    
def MSE_rootsolve(target_MSE,p,z,qt,T_guess,icemode=1):
    def obj_func(x):
        return target_MSE-MSE(T=x,z=z,qv=partition_mixedphase(p=p,T=x,qt=qt,icemode=icemode)['qv'],
                              ql=partition_mixedphase(p=p,T=x,qt=qt,icemode=icemode)['ql'],
                              qs=partition_mixedphase(p=p,T=x,qt=qt,icemode=icemode)['qs'])         
    T_solve = opt.fsolve(obj_func,T_guess)[0]
    water = partition_mixedphase(p,T_solve,qt=qt,icemode=icemode)
    qv = water['qv']
    ql = water['ql']
    qs = water['qs']
    if np.abs((MSE(T_solve,z,qv,ql,qs,icemode=icemode)-target_MSE)/target_MSE) > .001:
        print 'warning: MSE rootsolve failed'
        print target_MSE
        print p
        print z
        print qt
        print qv
        print ql
        print qs
        print T_solve
        return
    return T_solve
    
def adiabatic_profile(mode,T0,qt0,z0,p0,rho_e=0.,rho_z=0.,ztop=15000.,dz=10.,
                      fallout_factor=0.,icemode=1):
    
    # print "Arguments:"
    # print mode,T0,qt0,z0,p0,rho_e,rho_z,ztop,dz,fallout_factor,icemode
    # height
    z = np.arange(z0,ztop+dz,dz)
    
    # reference density (interpolated log-linearly to new z grid)
    if (mode=='parcel'):
        rho_e2 = np.exp(np.interp(z,rho_z,np.log(rho_e)))
    
    # moist static energy
    mse = np.empty(len(z))
    
    # pressure
    p    = np.zeros(len(z))
    logp = np.empty(len(z))
    
    # temperature, density, water mass fractions, gas constant, buoyancy
    T   = np.zeros(len(z))
    rho = np.zeros(len(z))
    b   = np.empty(len(z))
    qv  = np.empty(len(z))
    ql  = np.empty(len(z))
    qs  = np.empty(len(z))
    qt  = np.empty(len(z))
    Rm  = np.empty(len(z))
    
    water_init = partition_mixedphase(p0,T0,qt0,icemode=icemode)
    qv[0]      = water_init['qv']
    ql[0]      = water_init['ql']
    qs[0]      = water_init['qs']

    # initialize
    T[0]      = T0
    p[0]      = p0
    qt[0]     = qt0
    logp[0]   = np.log(p0)
    mse[0]    = MSE(T[0],z[0],qv[0],ql[0],qs[0],icemode=icemode)
    Rm[0]     = qv[0]*R_v + (1.-qv[0]-ql[0]-qs[0])*R_d
    rho[0]    = p[0]/(Rm[0]*T[0])
    if (mode == 'parcel'):
        b[0]  = gg*(rho_e2[0]/rho[0]-1.)
    elif (mode == 'env'):
        b[0] = 0.
        
    if (qv[0] >= qvstar(p[0],T[0],ql=ql[0],qs=qs[0],icemode=icemode)):
        unsat = False
    else:
        unsat = True
    
    for i in range(1,len(z)):
        # calculate pressure at next level
        if (mode == 'parcel'):
        	logp[i] = logp[i-1] - dz*gg*rho_e2[i-1]/p[i-1]
        	p[i] = np.exp(logp[i])
        elif (mode == 'env'):
            logp[i] = logp[i-1] - dz*gg/(Rm[i-1]*T[i-1])
            p[i] = np.exp(logp[i])

        ## ADDED FOR DEBUG
        if p[i] == 0:
            continue

        # calculate MSE at next level
        mse[i] = mse[i-1] - b[i-1]*dz
        
        if (unsat):
            cpm = qv[i-1]*c_pv + (1.-qv[i-1])*c_pd
            T[i] = (mse[i] - gg*z[i] - qv[i-1]*(E_0v + R_v*T_0) + cpm*T_0)/cpm
            
            if qvstar(p[i],T[i],icemode=icemode) <= qv[i-1]:
                #print 'Saturation at ', z[i],'m'
                unsat = False

        else:
            T[i] = MSE_rootsolve(mse[i],p[i],z[i],qt[i-1],T[i-1],icemode=icemode)
            

        water = partition_mixedphase(p[i],T[i],qt=qt[i-1],icemode=icemode)
        qv[i] = water['qv']
        ql[i] = water['ql']
        qs[i] = water['qs']
        
        if (fallout_factor > 0.):
                adjusted_state = fallout(fallout_factor,p[i],T[i],
                                         qv[(i-1):(i+1)],ql[(i-1):(i+1)],qs[(i-1):(i+1)],
                                         icemode=icemode,verbose=False)
                ql[i] = adjusted_state['ql_new']
                qs[i] = adjusted_state['qs_new']
                mse[i] = MSE(T[i],z[i],qv[i],ql[i],qs[i],icemode=icemode)
        
        qt[i] = qv[i] + ql[i] + qs[i]
        
        Rm[i] = qv[i]*R_v + (1.-qv[i]-ql[i]-qs[i])*R_d
        rho[i] = p[i]/(Rm[i]*T[i])
        if (mode=='parcel'):
            b[i] = gg*(rho_e2[i]/rho[i] - 1.)
        elif (mode=='env'):
            b[i] = 0.
        
    return_dict = {}
    return_dict['MSE'] = mse
    return_dict['T'] = T
    return_dict['b'] = b
    return_dict['rho'] = rho
    return_dict['qv'] = qv
    return_dict['ql'] = ql
    return_dict['qs'] = qs
    return_dict['qt'] = qt
    return_dict['p']  = p
    return_dict['z'] = z
    return return_dict
    
def adiabat_simple(T_base,qt,zbot,ztop,dz,p,z_les):
    z = np.arange(zbot,ztop+dz,dz)
    p = np.exp(np.interp(z,z_les,np.log(p))) 
    
    qv    =  np.zeros(len(z))      # kg/kg     specific humidity of plume
    qc    =  np.zeros(len(z))      # kg/kg     condensate (liquid only) mass fraction of plume
    qsat  =  np.zeros(len(z))      # kg/kg     saturation specific humidity of plume/env
    # Temperature
    T     =  np.zeros(len(z))      # K         
    # Density
    rho   =  np.zeros(len(z))      # kg/m^3 
    
    # Check if there is any liquid or solid water at plume base
    if qt > qvstar(p[0],T_base,liquid=True):
        print 'Warning: plume base supersaturated. Code not verified.'
        
    T[0]   = T_base
    qv[0]  = qt
    qc[0]  = 0.
        
    MSE = MSE_simple(z[0],T_base,qt)
    
    # Flag for LCL z index
    LCL = 0
    
    # Integrate plume upwards in height
    for i in range(len(z)):
        
        # Calculate plume, environment density
        rho[i] = p[i]/(R_d*T[i])
        
        # Calculate plume saturation specific humidity
        qsat[i]  = qvstar(p[i],T[i],liquid=True)
        
        if (((qt - qsat[i]) > 0) and (LCL == 0)):
                LCL = 1
            
        if i < (len(z)-1):
            # Calculate temperature
            # Use root finding algorithm fsolve with initial guess
            # Define objective function (redefined at each step) for root solver
            def rootsolve_MSE(target_MSE,z,p,T_upperlim,qv_upperlim):
                def obj_func(x):
                    return target_MSE-MSE_simple(z,x,qvstar(p,x,liquid=True))
                T_solve = opt.fsolve(obj_func,T_upperlim)[0]
                if T_solve >= T_upperlim:
                    print 'Warning: T not advanced'
                    print T_solve
                if qvstar(p,T_solve,liquid=True) >= qv_upperlim:
                    print 'Warning: invalid vapor mixing ratio'
                    print i, qvstar(p,T_solve,liquid=True), qv_upperlim
                return T_solve
 
            if LCL == 1:
                T[i+1] = rootsolve_MSE(MSE,z[i+1],p[i+1],T[i],qt)
            else:
                T[i+1] = (MSE - gg*z[i+1] - L_v*qv[i])/c_pd
        
        	# partition water
            qv[i+1] = np.min([qvstar(p[i+1],T[i+1],liquid=True),qt])
            qc[i+1] = qt - qv[i+1] 
        
    return_dict = {}
    return_dict['rho'] = rho
    return_dict['T'] = T
    return_dict['qv'] = qv
    return_dict['qc'] = qc
    return_dict['p'] = p
    return_dict['z'] = z
    return_dict['MSE'] = MSE
    return (return_dict)

