




###--- Modules ---###
from math import *
import numpy as np
# import sys


###--- Functions for rainfall statistics ---###
    
## Define inverse-logarithmic percentiles and bin edges
## Return: maximum percentile, percentiles (independent from values)
##     and centers_invlog (information about values)
## n_unit_min is the minimum number of values per bin, used to compute highest percentile
def defineInvLogPercentiles(values,subset=None,n_unit_min=4):
	# Reduce values to subset
	if subset is None:
		n_pts = values.size
		vals = values
	else:
		n_pts = subset.sum()
		sshape = subset.shape
		if len(sshape) == 2:
			ntime = values.shape[0]
			ind_sub = np.vstack([[subset]]*ntime)
		elif len(sshape) == 3:
			ind_sub = subset
		vals = values[ind_sub]
	n_max = int(log10(n_pts/n_unit_min))
	di = 0.1
	max_percentile = (1.-10**(-n_max))*100
	scale_invlog = np.arange(0,n_max+di,di)
	perc_invlog = np.subtract(np.ones(scale_invlog.size),np.power(10,-scale_invlog))*100
	centers_invlog = np.percentile(vals,perc_invlog)
	breaks_invlog = (centers_invlog[1:]+centers_invlog[:-1])/2
	binsize_invlog = np.diff(breaks_invlog)
	centers_invlog = centers_invlog[1:-1]
	perc_invlog = perc_invlog[1:-1]
	return max_percentile, perc_invlog, centers_invlog, breaks_invlog

## Define logarithmic bin edges accordingly to a set of values
## Return edges
def defineLogBins(values,subset=None,vmin=None):
	di = 0.1
	# Reduce values to subset
	if subset is None:
		vals = values
	else:
		sshape = subset.shape
		if len(sshape) == 2:
			ntime = values.shape[0]
			ind_sub = np.vstack([[subset]]*ntime)
		elif len(sshape) == 3:
			ind_sub = subset
		vals = values[ind_sub]
	# Compute vmin		
	if vmin is None:
		vals_nonzero = vals.copy()
		vals_nonzero[vals_nonzero==0] = None
		vmin = abs(np.nanmin(vals_nonzero))
	# Calculate log bins
	expmin = floor(log10(vmin))
	expmax = ceil(log10(vals.max()))
	exps = np.arange(expmin,expmax,di)
	breaks_log = np.power(10.,exps)
	centers_log = 0.5*(breaks_log[:-1]+breaks_log[1:])
	return centers_log, breaks_log

## Compute inverse-log percentiles, values and edges,
## and logarithmic values, edges and density
def getStatistics(values,subset=None,n_unit_min=4):
	# Reduce values to subset
	if subset is None:
		vals = values
	else:
		sshape = subset.shape
		if len(sshape) == 2:
			ntime = values.shape[0]
			ind_sub = np.vstack([[subset]]*ntime)
		elif len(sshape) == 3:
			ind_sub = subset
		vals = values[ind_sub]
	print "Sample size :", vals.size
	# Compute statistics
	maxQ, Q_IL, vals_c_IL, vals_e_IL = defineInvLogPercentiles(vals,n_unit_min=n_unit_min)    # For "Quantiles Inverse-Logarithmic" 
	vals_c_L, vals_e_L = defineLogBins(vals)    # For "edges" and "centers", "Logarithmic"
	vals_Hd_L, vals_e_L = np.histogram(vals,bins=vals_e_L,density=True)    # For "Histogram density"
	vals_Hd_L[vals_Hd_L==0] = None
	vals_Hd_IL, vals_e_IL = np.histogram(vals,bins=vals_e_IL,density=True)    # For "Histogram density"
	vals_Hd_IL[vals_Hd_IL==0] = None
	return maxQ, Q_IL, vals_c_IL, vals_e_IL, vals_c_L, vals_e_L, vals_Hd_L, vals_Hd_IL

## Compute 2D joint statistics on inverse-logarithmic axes
def getJointStatistics(values1,values2,subset=None,n_unit_min=4):
	# Reduce values to subset
	# Define ind_sub
	ntime = values1.shape[0]
	if subset is None:
		ind_sub = np.vstack([[subset]]*ntime)
	else:
		sshape = subset.shape
		if len(sshape) == 2:			
			ind_sub = np.vstack([[subset]]*ntime)
		elif len(sshape) == 3:
			ind_sub = subset
	if subset is None:
		if values1 is None or values2 is None:
			return
		vals1 = values1.ravel()
		vals2 = values2.ravel()
	else:
		vals1 = values1[ind_sub].ravel()
		vals2 = values2[ind_sub].ravel()
	print "Sample size :", vals1.size, vals2.size
	# Compute statistics
	maxQ1, Q1_IL, vals1_c_IL, vals1_e_IL = defineInvLogPercentiles(vals1,n_unit_min=n_unit_min)
	maxQ2, Q2_IL, vals2_c_IL, vals2_e_IL = defineInvLogPercentiles(vals2,n_unit_min=n_unit_min)
	vals_H2Dd_IL, vals1_e_IL, vals2_e_IL = np.histogram2d(x=vals1,y=vals2,bins=(vals1_e_IL,vals2_e_IL),normed=True)
	return maxQ1, maxQ2, Q1_IL, Q2_IL, vals1_e_IL, vals1_c_IL, vals2_e_IL, vals2_c_IL, vals_H2Dd_IL

## From values and percentile (in percents), find
def getIndicesOfExtremePercentile(values,percentile,subset=None,n_unit_min=4):
	if subset is None:
		vals = values
	else:
		# Reduce values to subset
		ntime = values.shape[0]
		sshape = subset.shape
		if len(sshape) == 2:
			ind_sub = np.vstack([[subset]]*ntime)
		elif len(sshape) == 3:
			ind_sub = subset
		vals = values[ind_sub]
	# Get percentiles in inverse-logarithmic fashion
	maxQ, Q_IL, vals_c_IL, vals_e_IL = defineInvLogPercentiles(vals,n_unit_min=n_unit_min)
	# Index of quantile as the index of closest value Q_IL to percentile
	dist_to_Q_IL = np.absolute(np.subtract(Q_IL,percentile*np.ones(Q_IL.shape)))
	mindist = dist_to_Q_IL.min()
	i_Q = np.argmax(dist_to_Q_IL == mindist)
	# Get limit values
	vmin = vals_e_IL[i_Q]
	if i_Q < len(Q_IL)-1:
		vmax = vals_c_IL[i_Q+1]
	else:
		vmax = vmin*1000.
	# Get index of locations
	ind_Q = np.logical_and(values >= vmin, values < vmax)
	if subset is not None:
		ind_Q = np.logical_and(ind_Q,ind_sub)
	# If Q is the last one, zero out the subset indices
	if i_Q == len(Q_IL)-1:
		ind_Q = np.zeros(ind_Q.shape,dtype=bool)
	return ind_Q

## Compute mean of values over points referred to by ind.
## Conserves vertical profile if necessary, if in second dimension.
def computeMeanAtTimeLatLonIndices(values,ind):
	if ind.sum() == 0:
		return None
	else:
		vshape = values.shape
		ishape = ind.shape
		sameShape = (vshape == ishape)
		sameShapePlusVertical = (len(vshape) == len(ishape)+1 and
								vshape[0] == ishape[0] and 
								vshape[-2:] == ishape[-2:])
		if sameShape:
			vmean = values[ind].mean()
		elif sameShapePlusVertical:
			nlev = vshape[1]
			ind_ok = np.vstack([[ind]]*nlev)
			values_lev_1st = values.swapaxes(0,1)
			vmean = values_lev_1st[ind_ok].reshape(nlev,ind.sum()).mean(axis=1)
		return vmean

## Compute time-horizontal mean over all points
def computeMean(values,subset=None):
	vshape = values.shape
	isTimeLatLon = (len(vshape) == 3)
	isTimeLevLatLon = (len(vshape) == 4)
	if isTimeLatLon:
		if subset is None:
			return values.mean()
		else:
			print type(subset)
			sshape = subset.shape
			if len(sshape) == 2:
				ntime = vshape[0]
				ind_subset = np.stack([subset]*ntime)
			elif len(sshape) == 3:
				ind_subset = subset
			return values[ind_subset].mean()
	elif isTimeLevLatLon:
		if subset is None:
			return np.mean(values,axis=(0,2,3))
		else:
			sshape = subset.shape
			nlev = vshape[1]
			if len(sshape) == 2:
				ntime = vshape[0]
				ind_subset = np.vstack([[np.vstack([[subset]]*ntime)]]*nlev)
				values_lev_1st = values.swapaxes(0,1)
				return values_lev_1st[ind_subset].reshape(nlev,ntime*subset.sum()).mean(axis=1)
			elif len(sshape) == 3:
				ind_subset = np.vstack([[subset]]*nlev)
				values_lev_1st = values.swapaxes(0,1)
				return values_lev_1st[ind_subset].reshape(nlev,subset.sum()).mean(axis=1)





