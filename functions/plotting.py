###############################################################################
#                                                                             #
#   Functions to plot maps, scatter plots, etc. compatible with swath and     #
#   gridded preprocessed data.                                                #
#                                                                             #
#                                                     Benjamin Fildier 2016   #
#                                                                             #
###############################################################################

from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as tk
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm
from scipy import stats
from math import *
from datetime import date, time, datetime
import string
import numpy as np
import numpy.ma as ma
import os, sys

## Add ../functions directory to system path
currentpath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,currentpath[:currentpath.rfind('/')+1]+'functions')
## Own libraries
# from correlationsAndScalograms import *

### Functions ###

## Flip longitudes between (0,360) and (-180,180) and vice-versa
def flipLongitudesToSymmetric(lons):
    lons[lons >= 180] = ma.add(lons[lons >= 180],\
                -360*ma.masked_array(np.ones(lons.shape), lons < 180)[lons >= 180])
    return lons
def flipLongitudesToPositive(lons):
    lons[lons < 0] = ma.add(lons[lons < 0],\
                360*ma.masked_array(np.ones(lons.shape), lons >= 0)[lons < 0])
    return lons

## Get string representing the (square) subdomain based on domain range
def getSubdomainString(lonlims,latlims):
	if latlims[0] < 0:
		subdomain_string = str(int(abs(latlims[0])))+'S'+str(int(lonlims[0]))+'E'
	else:
		subdomain_string = str(int(abs(latlims[0])))+'N'+str(int(lonlims[0]))+'E'
	dlon = str(int(lonlims[1]-lonlims[0]))
	dlat = str(int(latlims[1]-latlims[0]))
	subdomain_string = subdomain_string+'_'+dlon+'by'+dlat
	return subdomain_string

def initSubplot(lons_range,lats_range, varname, units, range_dates_times, ax, vmin=None, vmax=None, log=False, cmap=plt.cm.jet,ftsz=14):
	# create polar stereographic Basemap instance.
	m = Basemap(projection='cyl',lon_0=(lons_range[0]+lons_range[1])/2,lat_0=(lats_range[0]+lats_range[1])/2,lat_ts=0,\
				llcrnrlat=lats_range[0],urcrnrlat=lats_range[1],\
				llcrnrlon=lons_range[0],urcrnrlon=lons_range[1],\
				rsphere=6371200.,resolution='l',area_thresh=10000)
	m.drawcoastlines()
	fillconts = True
	if fillconts:
		m.fillcontinents()
	# derive spacing for meridians and parallels
	dlon = lons_range[1]-lons_range[0]
	if dlon < 50:
		p = log10(dlon)
		q = floor(p)
		dl = 10**q
		if p-q < 0.2:
			dl = dl/2.
		if p-q > 0.8:
			dl = dl*2.
	else:
		dl = 30
	## draw parallels
	parallels = np.arange(-90.,90.,dl)
	m.drawparallels(parallels,labels=[1,0,0,0],fontsize=ftsz-4)
	## draw meridians
	# meridians = np.arange(10*(lons_range[0]//10),10*(lons_range[1]//10+1),dl)
	meridians = np.arange(0,360,dl)
	m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=ftsz-4)
	if range_dates_times != None:
		#     if range_dates_times['start'].date() == range_dates_times['end'].date():
		#         times = range_dates_times['start'].isoformat(' ')+' to '+range_dates_times['end'].time().isoformat()
		#     else:
		#         times = range_dates_times['start'].isoformat(' ')+' to \n'+range_dates_times['end'].isoformat(' ')
		times = range_dates_times['equator'].isoformat(' ')
	else:
		times = ''
	## Define title
	if type(units) == type(None) and times == '':
		title = varname
	elif times == '':
		title = varname+' ('+units+')'
	elif type(units) == type(None):
		title = varname+'\n'+times
	else:
		title = varname+' ('+units+')\n'+times
	# ax.set_title(title,fontsize=ftsz)
	plt.title(title,fontsize=ftsz)
	return m

def addToSubplot(m,lons, lats, values, vmin=None, vmax=None, log=False, cmap=plt.cm.jet):
	if log:
		m.pcolormesh(lons,lats,values, shading='flat',cmap=cmap,latlon=True,norm=LogNorm(vmin=vmin,vmax=vmax))
	else:
		m.pcolormesh(lons,lats,values, shading='flat',cmap=cmap,latlon=True,vmin=vmin,vmax=vmax)

def drawSubplot(lons_range,lats_range, lons, lats, values, varname, units, range_dates_times, ax, vmin=None, vmax=None, log=False, cmap=plt.cm.jet,ftsz=14,cbar=True):
	## Initialize subplot
	m = initSubplot(lons_range,lats_range, varname, units, range_dates_times, ax, vmin, vmax, log, cmap,ftsz)
	## Draw data
	addToSubplot(m,lons, lats, values, vmin, vmax, log, cmap)
	## Add colorbar
	if cbar:
		cb = m.colorbar(location='right',size="5%",pad="2%")
		cb.ax.tick_params(labelsize=14)
	
def scatterDensitySubplot(x,y,xlab,ylab,xunits='',yunits='',xlim=(None,None),\
	ylim=(None,None),bins=None,title=None,cmap=plt.cm.Reds,ftsz=18):
	valid_points = commonPoints(x,y)
	plt.hexbin(x[valid_points],y[valid_points],bins=bins,cmap=cmap)
	plt.xlabel(xlab+" ("+xunits+")",fontsize=ftsz)
	plt.ylabel(ylab+" ("+yunits+")",fontsize=ftsz)
	plt.xlim(xlim)
	plt.ylim(ylim)
	cb = plt.colorbar()
	if title==None:
		if bins=='log':
			plt.title('Number of points (log10(N))',fontsize=ftsz)
			cb.set_label('Number of points '+bins+'(N)',fontsize=ftsz-6)
		else:
			plt.title('Number of points (N)',fontsize=ftsz)
			cb.set_label('Number of points (N)',fontsize=ftsz-6)
	else:
		plt.title(title,fontsize=ftsz)

## Draw scatter plot with given labels and title
def scatterSubplot(x,y,xlab,ylab,xunits='',yunits='',title='',color=None):
	plt.scatter(x,y,c=color)
	plt.xlabel(xlab+" ("+xunits+")",fontsize=16)
	plt.ylabel(ylab+" ("+yunits+")",fontsize=16)
	plt.title(title)
## Draw histogram of the data
def histogramSubplot(x,bins=10,range=None,xlab='',xunits='',title=''):
	plt.hist(x,bins=bins,range=range)
	plt.xlabel(xlab+" ("+xunits+")",fontsize=16)
	plt.ylabel('Probability Density',fontsize=16)
	plt.title(title)

### Error model ###

## Plot the sigma error as a function of the fraction of missing data and the coarsening level
## 'modelNP' is the sigma error model function that takes an array of numpy arrays as the first argument
def derivationSigmaErrorPlot(error_sigma,f_values,levels,modelNP,coefs,title='',r_index=0):
	fig = plt.figure(figsize=(14,6))
	fig.suptitle(title,fontsize=20)
	gs = mpl.gridspec.GridSpec(1, 2, width_ratios=[6, 7]) 
	nf = len(f_values)
	nl = len(levels)
	err_range = error_sigma[:,:,r_index].min(),error_sigma[:,:,r_index].max()
	derr = err_range[1]-err_range[0]
	ax = fig.add_subplot(gs[0])
	ax.tick_params(axis='both', which='major', labelsize=13)
	legends = ['Level 2: 4x4 Points','Level 3: 8x8 Points','Level 4: 16x16 Points',\
          'Level 5: 32x32 Points']
	h = []
	colors = ('r','k','b','g')
	for i in range(nl):
		h.append(plt.scatter(f_values,error_sigma[i,:,r_index],color=colors[i]))
		plt.plot(f_values,modelNP(np.array([[levels[i]]*nf,f_values]),*coefs),\
				dashes=[2, 2, 2, 2],color=colors[i])
	plt.xlabel('Fraction Missing '+r'$f$',fontsize=20)
	plt.ylabel('Sigma Error '+r'$\sigma_\varepsilon$',fontsize=20)
	plt.ylim(err_range[0]-derr/10.,err_range[1]+derr/10.)
	plt.legend(tuple(h),tuple(legends),scatterpoints=1,loc='upper left',ncol=1,fontsize=12)
	ax = fig.add_subplot(gs[1])
	ax.tick_params(axis='both', which='major', labelsize=13)
	sm = plt.cm.ScalarMappable(cmap=plt.cm.brg,norm=mpl.colors.Normalize(vmin=f_values.min(),vmax=f_values.max()))
	for i in range(nf):
		col = sm.to_rgba((f_values[i]))
		if nl == 4:
			col = col[:3]
		plt.scatter(levels,error_sigma[:,i,r_index],color=col,edgecolors='none')
		plt.plot(levels,modelNP(np.array([levels,[f_values[i]]*nl]),*coefs),\
				dashes=[2, 2, 2, 2],color=col)
	plt.xlabel('Coarsening Level, '+r'$l = \frac{1}{2} \log_2(n_l)$',fontsize=20)
	plt.ylabel('Sigma Error '+r'$\sigma_\varepsilon$',fontsize=20)
	plt.ylim(err_range[0]-derr/10.,err_range[1]+derr/10.)
	sm._A = []    	# Mysterious command necessary to use it as argument in colorbar 
	 			  	# Otherwise get TypeError: You must first set_array for mappable)
	cb = plt.colorbar(sm)
	cb.set_label('Fraction Missing '+r'$f$',fontsize=16)

## Plot the Pearson coefficient as a function of the scale
def pearsonRVsScaleSubplot(pearson_coefs,levels,resolution,error,xlab,ylab,title=''):
	x = np.power(2,levels)*resolution
	dx = x.max()-x.min()
	y1 = np.add(pearson_coefs,error)
	y2 = np.subtract(pearson_coefs,error)
	plt.semilogx(x,pearson_coefs,marker='o',c='g')
	plt.xlim((x.min()-dx/10.,x.max()+dx/10.))
	plt.ylim((0,1))
	plt.fill_between(x, y1, y2, facecolor='green', interpolate=True, alpha=0.5)
	plt.xlabel(xlab,fontsize=14)
	plt.ylabel(ylab,fontsize=14)
	plt.title(title,fontsize=16)

## Draw subplot of values (pearson coef or absolute errors on it) on a 2D plot as a 
## function of the coaserning level (x axis) and the maximum allowed fraction of 
## missing data f
def valuesVsScaleAndFmaxSubplot(values,levels,resolution,fmaxs,n_values,title='',cmap=plt.cm.jet,ax=None,nmin=0,vmin=0,vmax=None):
	nf = len(fmaxs)
	nl = len(levels)
	fmax_grid, levels_grid = np.meshgrid(fmaxs,levels)
	new_values = ma.masked_array(values,n_values<nmin)
	values_to_plot = np.swapaxes(np.fliplr(new_values),0,1)
	im = ax.imshow(values_to_plot,cmap=cmap,aspect='auto',interpolation='none',vmin=vmin,vmax=vmax)
	ax.tick_params(axis='both', direction='out', labelsize=14)
	ax.set_yticks(np.arange(nf))
	ax.set_yticklabels(fmaxs.astype('str')[::-1])
	ax.set_xticks(np.arange(nl))
	xticklabs = map(lambda x:'%1g'%x,resolution*np.power(2.,levels))
	ax.set_xticklabels(xticklabs)
	plt.xlabel('Scale (km)',fontsize=14)
	plt.ylabel('Maximum f',fontsize=14)
	plt.colorbar(im)
	plt.title(title,fontsize=16)


## Draw subplot of structure function as a function of scale for various variables
def structureFunctionsSubplot(strfns,dk_values,resolution,varnames,errors=None,xlab='',\
	ylab='',title='',l_min=0,order=1,ax=None,ftsz=18):
	x = dk_values*resolution
	# Define colors
	nvar = len(strfns)
	sm = plt.cm.ScalarMappable(cmap=plt.cm.gnuplot2,norm=mpl.colors.Normalize(vmin=0,vmax=nvar))
	xmin = x.mean()
	xmax = x.mean()
	ymin = 1.
	ymax = 1.
	for i in range(0,nvar):
		il_ref = np.where(dk_values == 2**l_min)
		factor = 1./np.nanmin(strfns[i][il_ref])
		col = sm.to_rgba(i)
		valid = np.logical_and(strfns[i]!=0,dk_values >= 2.**l_min)
		if varnames[i] == 'Wind Speed':
			valid = np.logical_and(valid,dk_values >= 3)
		xmin = min(xmin,x[valid].min())
		xmax = max(xmax,x[valid].max())
		ymin = min(ymin,np.nanmin(np.array(strfns[i][valid]))*factor)
		ymax = max(ymax,np.nanmax(np.array(strfns[i][valid]))*factor)
		plt.loglog(x[valid],strfns[i][valid]*factor,marker='o',color=col,label=varnames[i])
		if type(errors) != type(None):
			y1 = np.array(strfns[i],dtype=np.float) + np.array(errors[i],dtype=np.float)
			valid = np.logical_and(valid,np.isnan(y1))
			y2 = np.array(strfns[i],dtype=np.float) - np.array(errors[i],dtype=np.float)
			y1[y1 <= 0.] = 1e-5
			y2[y2 <= 0.] = 1e-5
			plt.fill_between(x[valid], y1[valid]*factor, y2[valid]*factor, facecolor=col, interpolate=True, alpha=0.3)
	ax.set_xlim((xmin/1.2,xmax*1.2))
	ax.set_ylim((ymin/1.1,ymax*1.2))
	plt.legend(scatterpoints=1,loc='upper left',ncol=2,fontsize=ftsz-4)
	if xlab == '':
		plt.xlabel('Scale (km)',fontsize=ftsz-2)
	else:
		plt.xlabel(xlab,fontsize=ftsz-2)
	if ylab == '':
		plt.ylabel('Structure Function',fontsize=ftsz-2)
	else:
		plt.ylabel(ylab,fontsize=ftsz-2)
	if title == '':
		plt.title('Order q=%g'%order,fontsize=ftsz)
	else:
		plt.title(title,fontsize=ftsz)

## Draw subplot of values (pearson coef or absolute errors on it) as a function
## of the coaserning level (x axis) and overlay the range of errors
def pearsonRVsScaleAllFmaxSubplot(pearson_coefs,errors,n_values,levels,resolution,fmaxs,\
	ylab='',title='',cmap=plt.cm.jet,ax=None,nmin=0,rmin=0,rmax=1,vmin=0,vmax=None,ftsz=16,cbar=True):
	x = np.power(2,levels)*resolution
	dx = x.max()-x.min()
	nfmax = len(fmaxs)
	sm = plt.cm.ScalarMappable(cmap=cmap,norm=mpl.colors.Normalize(vmin=vmin,vmax=vmax))
	for ifmax in range(nfmax)[::-1]:
		fmax = fmaxs[ifmax]
		col = sm.to_rgba(fmax)
		mask_n = n_values[:,ifmax] < nmin
		y = ma.masked_array(pearson_coefs[:,ifmax],mask_n)
		plt.semilogx(x,y,marker='.',c=col,ls='--',linewidth=3.)
		dy = ma.masked_array(errors[:,ifmax],n_values[:,ifmax] < nmin)
		y1 = np.array(y,dtype=float) + np.array(dy,dtype=float)
		y1 = ma.masked_array(y1,mask_n)
		y2 = np.array(y,dtype=float) - np.array(dy,dtype=float)
		y2 = ma.masked_array(y2,mask_n)
		plt.fill_between(x, y1, y2, facecolor=col, interpolate=True, alpha=0.3, linewidth=0.)
	xlim = (x.min()/1.2,x.max()*1.2)
	ax.add_line(Line2D((x.min()/10.,x.max()*10.),[0.,0.],c='k',linestyle='--',linewidth=1.))
	ax.set_xlim(xlim)
	plt.ylim((rmin,rmax))
	if ylab == '':
		ylab = 'Pearson coefficient R'
	plt.xlabel('Scale (km)',fontsize=ftsz-2)
	plt.ylabel(ylab,fontsize=ftsz-2)
	plt.title(title,fontsize=ftsz+2)
	sm._A = []    	# Mysterious command necessary to use it as argument in colorbar 
	 			  	# Otherwise get TypeError: You must first set_array for mappable)
	if cbar:
		cb = plt.colorbar(sm)
		cb.set_label('Maximum value allowed for $f$',fontsize=ftsz-8)









