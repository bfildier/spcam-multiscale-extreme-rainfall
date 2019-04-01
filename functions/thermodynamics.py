##################
# Thermodynamics #
##################

## Modules
from math import *
import numpy as np
import numpy.ma as ma

## Own modules
# from physicalConstants import *
from thermo_constants import *
g = gg

###--- Some thermodynamic functions ---###

## Compute the air density profile based on profiles for T, p and q
## Use ideal gas law and virtual effect in linearized form
def airDensity(temp,pres,shum):
    rho_dry = np.divide(pres,R_d*temp)
    return np.divide(rho_dry,(1+eps*shum))

# Compute the saturation vapor pressure from expressions by Buck (1996)
def saturationVaporPressure(temp):
    T_0 = 273.15
    whereAreNans = np.isnan(temp)
    temp_wo_Nans = temp.copy()
    temp_wo_Nans[whereAreNans] = 0.
    # Initialize
    e_sat = np.zeros(temp.shape)
    e_sat[whereAreNans] = np.nan
    # T > 0C
    overliquid = np.array((temp_wo_Nans > T_0),dtype=bool)
    e_sat_overliquid = 0.61121*np.exp(np.multiply(18.678-(temp-T_0)/234.5,
                                                  np.divide((temp-T_0),257.14+(temp-T_0))))
    e_sat[overliquid] = e_sat_overliquid[overliquid]
    # T < 0C 
    overice = np.array((temp_wo_Nans < T_0),dtype=bool)
    e_sat_overice = 0.61115*np.exp(np.multiply(23.036-(temp-T_0)/333.7,
                                               np.divide((temp-T_0),279.82+(temp-T_0))))
    e_sat[overice] = e_sat_overice[overice]
    
    return e_sat*1000       # in Pa

# Compute the saturation specific humidity based on the expressions by Buck
def saturationSpecificHumidity(temp,pres):
    e_sat = saturationVaporPressure(temp)
    return np.divide(e_sat/R_v,pres/R_d)


## Compute z(p) using hydrostatic approximation
## rho is defined on p_vals axis; p_vals contains values of p_ref and p
def zOfP(p_ref,p,p_vals,rho,z_ref=0.):
    """" Calculates height of given pressure as the following:
    
    .. math::
        z(P) = \int_{p_{ref}}^{p} \\frac{dp}{\\rho(p)g}

    """
    i_ref = np.argmax(p_vals==p_ref)
    i_top = np.argmax(p_vals==p)

    if i_ref >= i_top:
        i_ref += 1
        z_values = z_ref + zOfP_all(p_vals[i_top:i_ref],rho[i_top:i_ref])
    else:
        i_top += 1
        i_top,i_ref = (i_ref,i_top)
        z_values = z_ref + zOfP_all(np.flipud(p_vals[i_top:i_ref]),
            np.flipud(rho[i_top:i_ref]))

    return z_values[0]


## Compute whole profile z(p) using hydrostatic approximation
def zOfP_all(p_vals,rho,z_ref=0.):
    if (len(p_vals) == 0 and len(rho) == 0):
        z_values = np.array([z_ref])
    else:
        i_ref = np.argmax(p_vals==p_vals.max())
        i_top = np.argmax(p_vals==p_vals.min())
        coef = 1
        if i_top > i_ref: 
            i_top,i_ref = (i_ref,i_top)
            coef = -1
        # Compute profile
        z = z_ref
        rho_mid = 0.5*(rho[i_top:i_ref]+rho[(i_top+1):(i_ref+1)])
        dp = np.diff(p_vals[i_top:(i_ref+1)])
        z_values = coef*np.divide(dp,g*rho_mid)
        if coef == 1:
            z_values = np.flipud(np.cumsum(np.append([z],np.flipud(z_values))))
        else:
            z_values = np.flipud(np.cumsum(np.append([z],z_values)))
    return z_values


## Compute p(z)
## rho is defined on z_vals axis; z_vals contains values of z_ref and z
def pOfZ(z_ref,z,z_vals,rho,p_ref):
    i_ref = np.argmax(z_vals==z_ref)
    i_top = np.argmax(z_vals==z)
    if i_ref <= i_top:    # if heights are in increasing order
        i_top += 1
        p_values = pOfZ_all(z_vals[i_ref:i_top],rho[i_ref:i_top],p_ref)
    else:    # if heights are in decreasing order
        i_ref += 1
        i_top,i_ref = (i_ref,i_top)
        p_values = pOfZ_all(z_vals[i_ref:i_top],rho[i_ref:i_top],p_ref)
    return p_values[0]


## Compute whole profile p(z) using hydrostatic approximation
def pOfZ_all(z_vals,rho,p_ref):
    if (len(z_vals) == 0 and len(rho) == 0):
        p_values = np.array([p_ref])
    else:
        i_ref = np.argmax(z_vals==z_vals.min())
        i_top = np.argmax(z_vals==z_vals.max())
        coef = 1
        if i_ref > i_top:
            i_top,i_ref = (i_ref,i_top)
            coef = -1
        # p = p_ref
        rho_mid = 0.5*(rho[i_ref:i_top]+rho[(i_ref+1):(i_top+1)])
        dz = np.diff(z_vals[i_ref:(i_top+1)])
        p_values = -coef*np.multiply(rho_mid,g*dz)
        if coef == 1:
            p_values = np.flipud(np.cumsum(np.append([p_ref],p_values)))
        else:
            p_values = np.flipud(np.cumsum(np.append([p_ref],np.flipud(p_values))))
    return p_values


###--- older functions ---###

# def saturatedVaporPressure(T): # in K
#    #  if isinstance(T,float):
#     #   return 6.1094*exp(17.625*(T-273.15)/(T-273.15+243.04))
#     # elif type(T) in [np.array,np.ma.core.MaskedArray]:
#     dims = T.shape
#     return 6.1094*np.exp(17.625*np.divide(np.add(T,-273.15*np.ones(dims)),\
#         np.add(T,(-273.15+243.04)*np.ones(T.shape))))

# ## Compute relative humidity from temperature and water vapor mixing ratio data
# def relativeHumidity(temp,wvmr,p=1000):   # p in hPa
#     if type(temp) != type(wvmr):
#         print 'The arguments should be both floats or both arrays of same dimensions.'
#     if type(temp) == float:
#         hspc = wvmr/1000/(wvmr/1000+1)
#         return p*hspc/((eps + (1-eps)*hspc)*saturatedVaporPressure(temp))
#     elif type(temp) == np.ndarray:
#         dims = temp.shape
#         hspc = np.divide(wvmr[:]/1000,np.add(np.ones(dims),wvmr[:]/1000))
#         sat_wv_vp = 6.1094*np.exp(17.625*np.divide(np.add(temp[:],-273.15*np.ones(dims)),\
#             np.add(temp[:],(-273.15+243.04)*np.ones(dims)))) # hPa
#         wv_vp = np.divide(p*hspc,np.add(eps*np.ones(dims),(1-eps)*hspc))
#         return np.divide(wv_vp,sat_wv_vp)  
#     elif type(temp) == np.ma.core.MaskedArray:
#         mask_common = np.logical_or(temp.mask,wvmr.mask)
#         dims = temp.shape
#         hspc = np.divide(wvmr[:]/1000,np.add(np.ones(dims),wvmr[:]/1000))
#         sat_wv_vp = 6.1094*np.exp(17.625*np.divide(np.add(temp[:],-273.15*np.ones(dims)),\
#             np.add(temp[:],(-273.15+243.04)*np.ones(dims)))) # hPa
#         wv_vp = np.divide(p*hspc,np.add(eps*np.ones(dims),(1-eps)*hspc))
# #         print "WV mixing ratio:", wvmr.min(), wvmr.max()
# #         print "WV vapor pressure:", wv_vp.min(), wv_vp.max()
# #         print "WV vapor pressure at saturation", sat_wv_vp.min(), sat_wv_vp.max()
#         hrel = np.divide(wv_vp,sat_wv_vp)  
#         return ma.masked_array(hrel,mask_common)
#     else:
#         print 'The arguments should be both floats or both arrays of same dimensions.'