
    



###--- Modules ---###
from netCDF4 import Dataset
from scipy import interpolate

## Own modules
from statistics import *
from thermodynamics import *
from thermo_constants import *
from thermo_funcs import *
g = gg

###--- Functions for scaling ---###

# Compute qvstar based on mean temperature profile
def qvstarProfile(pres,temp):

    """Compute profile of saturation specific humidity from pressure levels and
    corresponding temperature profile"""

    n = len(temp)
    qvs = np.zeros(n)
    for i in range(n):
        qvs[i] = qvstar(pres[i],temp[i],icemode=0.)
    return qvs

## profiles is a dictionary containing at least 'T' and 'p'
## At call, provide either p_vals or z_vals
def interpTempProfile(profiles,coord_vals,coord_name='z'):

    """profiles is a dictionary containing at least 'T' and 'p'
    At call, provide either p_vals or z_vals
    Function computeScaling using mode adiab_parcel would probably work
    only when interpolating with pressure coordinate"""

    if coord_name in ('z','p') and coord_name in profiles.keys():
        coord_old = profiles[coord_name]
        new_profiles = {coord_name:coord_vals}
        for field in ('T',):
            # print coord_old, profiles[field]
            new_profiles[field] = interpolate.interp1d(coord_old,profiles[field],kind='cubic')(coord_vals)
    else:
        print"'coord_name' must be 'z' or 'p' and must be present in 'profiles.keys()'"
        sys.exit(-1)
    return new_profiles

## Get slice for troposphere from mean temperature profile
def getTroposphereSlice(temp):

    """Compute a slice to crop the vertical profiles to relevant values.
    Here crop between the first level and the tropopoause defined as the height of minimum 'temp'."""

    if len(temp.shape) == 1:
        i_trop = np.argmax(temp == temp.min())
        i_BL = np.argmax(temp == temp.max())
    else:
        temp_m = computeMean(temp)
        i_trop = np.argmax(temp_m == temp_m.min())
        i_BL = np.argmax(temp_m == temp_m.max())
    i_min,i_max = i_trop,i_BL
    if i_min > i_max:
        i_min,i_max = i_max,i_min
    return slice(i_min,i_max+1)

## Crop all profiles entered as a tuple
def cropProfiles(profile_tuple,crop_slice):

    """Crop vertical profiles between boundary layer and tropopause if crop_slice
    is computed using function getTroposphereSlice."""

    cropped = ()
    for profile in profile_tuple:
        cropped += (profile[crop_slice],)
    return cropped


## Compute profiles to use for scaling based on mean profiles at quantile
def profilesForScalingFromMeanProfilesAtQ(pres,temp,omega,rho=None,w=None,qvs=None):

    # Compute qvstar profile according to reference temperature profile chosen
    if qvs is None:
        qvs = qvstarProfile(pres,temp)

    # Crop all profiles to tropospheric values
    slice_trop = getTroposphereSlice(temp)
    (pres_c,temp_c,qvs_c,omega_c) = cropProfiles((pres,temp,qvs,omega),slice_trop)
    if rho is not None:
        (rho_c,) = cropProfiles((rho,),slice_trop)
        if w is not None:
            (w_c,) = cropProfiles((w,),slice_trop)
            return pres_c, temp_c, qvs_c, omega_c, rho_c, w_c
        return pres_c, temp_c, qvs_c, omega_c, rho_c
    return pres_c, temp_c, qvs_c, omega_c

## Compute profiles to use for scaling based on mean profiles at quantile
# Qvstar computed with Buck's approximation
def profilesForScalingFromMeanProfilesAtQFromBuck(pres,temp,omega,rho=None,w=None):

    # Compute qvstar profile according to reference temperature profile chosen
    qvs = saturationSpecificHumidity(temp,pres)

    # Crop all profiles to tropospheric values
    slice_trop = getTroposphereSlice(temp)
    (pres_c,temp_c,qvs_c,omega_c) = cropProfiles((pres,temp,qvs,omega),slice_trop)
    if rho is not None:
        (rho_c,) = cropProfiles((rho,),slice_trop)
        if w is not None:
            (w_c,) = cropProfiles((w,),slice_trop)
            return pres_c, temp_c, qvs_c, omega_c, rho_c, w_c
        return pres_c, temp_c, qvs_c, omega_c, rho_c
    return pres_c, temp_c, qvs_c, omega_c

## Compute profiles to use for scaling based on mean profiles at quantile
def profilesForScaling2FromMeanProfilesAtQ(pres,temp,qvs,omega):

    # Crop all profiles to tropospheric values
    slice_trop = getTroposphereSlice(temp)
    (pres_c,temp_c,qvs_c,omega_c) = cropProfiles((pres,temp,qvs,omega),slice_trop)
    
    return pres_c,temp_c,qvs_c, omega_c


## Compute profiles later used for the scaling of extreme precipitation
def computeProfilesForScalingFromAllValues(mode,Q,pr,pres_surf,temp,spec_hum,omega,computeP,subset=None,n_unit_min=4):

    """Compute the vertical profiles required to compute O'Gorman and Schneider's scaling
    for precipitation extremes.
    'mode' can be 'global_mean', 'quantile_mean' or 'adiab_parcel' and refers to the 
    way of computing the reference temperature profile.
    All variable values are given in time,lat,lon coordinates."""
    
    # Get indices for a given quantile
    i_q = getIndicesOfExtremePercentile(pr,Q,subset=subset,n_unit_min=n_unit_min)
    # Define mean variables and profiles for this precipitation quantile
    ps_q = computeMeanAtTimeLatLonIndices(pres_surf,i_q)
    p_q = computeP(ps_q)
    t_q = computeMeanAtTimeLatLonIndices(temp,i_q)
    q_q = computeMeanAtTimeLatLonIndices(spec_hum,i_q)
    rho_q = airDensity(t_q,p_q,q_q)
    omega_q = computeMeanAtTimeLatLonIndices(omega,i_q)
    # pr_q = computeMeanAtTimeLatLonIndices(pr,i_q)

    # Choose and compute reference temperature profile
    temp_type = mode
    # temp_type = 'global_mean'
    if temp_type == 'global_mean':    # Global Mean T
        temp_ref = computeMean(temp)
    elif temp_type == 'quantile_mean':    # Mean T at locations of extreme quantile
        temp_ref = t_q
    elif temp_type == 'adiab_parcel':    # Use Jake's code
        i_BL = -1; i_top = 0
        z = zOfP_all(p_vals=p_q,rho=rho_q,z_ref=0.)
        dz = (z[i_top]-z[i_BL])/999.
        qt0 = qvstar(p_q[i_BL],t_q[i_BL])
        adiab_profiles = adiabatic_profile(mode='parcel',T0=t_q[i_BL],qt0=qt0,z0=z[i_BL],p0=p_q[i_BL],
                                       dz=dz,ztop=z[i_top],fallout_factor=1.,
                                       rho_e=rho_q,rho_z=z,icemode=1)
        adiab_profile_p = interpTempProfile(profiles=adiab_profiles,coord_vals=p_q,coord_name='p')
        temp_ref = adiab_profile_p['T']

    # Compute qvstar profile according to reference temperature profile chosen
    qvs_q = qvstarProfile(p_q,temp_ref)

    # Crop all profiles to tropospheric values
    slice_trop = getTroposphereSlice(temp)
    (p_c,qvs_c,omega_c) = cropProfiles((p_q,qvs_q,omega_q),slice_trop)

    return omega_c, qvs_c, p_c

## Compute vertical integral on pressure coordinates
def verticalPressureIntegral(pres,values=None):
    ## If values not given, compute weight of vertical column, normalizing factor
    if values is None:
        values = np.ones(pres.shape)
    ## Compute integral
    dp = np.diff(pres)
    val_mids = 0.5*(np.add(values[:-1],values[1:]))
    out = np.nansum(np.multiply(val_mids,dp)/g)
    return out

## Compute vertical integral of a product of two quantities on pressure coordinates
def verticalPressureIntegralProduct(pres,values1,values2=None,valmids2=None):
    dp = np.diff(pres)
    val_mids1 = 0.5*(np.add(values1[:-1],values1[1:]))
    if values2 is not None:
        val_mids2 = 0.5*(np.add(values2[:-1],values2[1:]))
    elif valmids2 is not None:
        val_mids2 = valmids2
    out = np.nansum(np.multiply(np.multiply(val_mids1,val_mids2),dp)/g)
    return out

## Compute scaling for extreme percentiles of precipitation based on vertical profiles
def computeScaling(omega,qvs=None,dqvsdp=None,pres=None):

    """From vertical profiles of omega and qvs, compute sum(omega*dqvs/g)"""

    if pres is None and qvs is not None:    # qvs must have the same dimension as omega

        # sign = np.sign(np.diff(qvs))[0]
        sign = 1
        dqvs = np.diff(qvs)
        omega_mid = 0.5*(np.add(omega[:-1],omega[1:]))
        # Compute integral + convert from Pa/(m/s) to mm/day
        pr_est = np.nansum(-sign*np.multiply(omega_mid,dqvs)/g)*86400
        
    elif pres is not None and dqvsdp is not None:    # qvs.size must equal omega.size-1

        dp = np.diff(pres)
        omega_mid = 0.5*(np.add(omega[:-1],omega[1:]))
        pr_est = np.nansum(-np.multiply(np.multiply(omega_mid,dqvsdp),
                                dp)/g)*86400

    return pr_est
