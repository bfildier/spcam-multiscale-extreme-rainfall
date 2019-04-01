################################################################################
#                                                                              #
#    Benjamin Fildier, July 2017.                                              #
#                                                                              #
#    From CRM hourly values of T and GCM hourly values of Ps, compute vertical #
#    profile of saturation specific humidity in each GCM grid cell, name it    #
#    CRM-QVSATENV (because use the environmental temperature profile CRM-T     #
#    to compute it in each CRM column before averaging).                       #
#                                                                              #
################################################################################


#### THIS SCRIPT ABORTS WHEN CALCULATING QVSTAR IN EACH CRM COLUMN BECAUSE OF
#### THE SIZE OF THE DATA TO OPERATE

###--- Modules ---###
import sys,os,glob
from netCDF4 import Dataset
import numpy as np
import socket
from argparse import ArgumentParser
import datetime as dt
from datetime import datetime,timedelta

## Own functions
currentpath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,currentpath[:currentpath.rfind('/')+1]+'functions')

from importingData import *
from thermodynamics import saturationSpecificHumidity,saturationVaporPressure


##--- Functions ---###

def getInputDirectories(dataroot,compset,experiment):

	hostname = socket.gethostname()
	case = "bf_%s_%s"%(compset,experiment)
	if hostname == "jollyjumper":
		inputdir = os.path.join(dataroot,"simulations",case)
		inputdir_processed_day = os.path.join(dataroot,'preprocessed',case,'day')
		inputdir_processed_1hr = os.path.join(dataroot,'preprocessed',case,'1hr')
	elif "edison" in hostname:
		inputdir = os.path.join(dataroot,'archive',case,"atm/hist")
		inputdir_processed_day = os.path.join(os.path.dirname(currentpath),'preprocessed',
			case,'day')
		inputdir_processed_1hr = os.path.join(os.path.dirname(currentpath),'preprocessed',
			case,'1hr')
	inputdir_results = os.path.join(os.path.dirname(currentpath),'results')
	inputdir_fx = os.path.join(dataroot,'preprocessed/allExperiments/fx')

	return inputdir, inputdir_processed_day, inputdir_processed_1hr,\
		inputdir_results, inputdir_fx

def parseArguments():

	parser = ArgumentParser(description="Extract GCM points corresponding to "+\
		"a given quantile of precipitation.")
	parser.add_argument('-f','--fraction',required=True,
		help="Threshold fraction alpha of rainfall amount to define rainfall event (%)."+\
			"If the integer entered is negative, all points are taken")
	parser.add_argument('-d','--dataroot',required=True,
		help="Input directory root.")
	parser.add_argument('-e','--experiment',default="piControl",
		choices=['piControl','abrupt4xCO2'],
		help="Physical experiment or scenario.")
	parser.add_argument('-c','--compset',choices=['FSPCAMm_AMIP','FAMIPC5'],
		default="FSPCAMm_AMIP",help="CESM compset used.")
	parser.add_argument('-dt','--datetimes',nargs=2,required=True,
		help="Boundary datetime values in YYYY-MM-DD-SSSSS format.")
	args = parser.parse_args()

	# Get time boundaries
	dt_bnds = []
	for tstamp in args.datetimes:
		dt_info = [int(s) for s in tstamp.split('-') if s.isdigit()]
		dt_bnds.append(datetime(dt_info[0],dt_info[1],dt_info[2])+timedelta(seconds=dt_info[3]))
	# Change date format to CMIP5 standard
	datetime1 = dt_bnds[0].isoformat().replace('-','').replace('T','').replace(':','')[:-2]
	datetime2 = dt_bnds[1].isoformat().replace('-','').replace('T','').replace(':','')[:-2]

	# Get case name
	case = "bf_%s_%s"%(args.compset,args.experiment)

	# Get input directories
	inputdir, inputdir_processed_day, inputdir_processed_1hr, inputdir_results, \
		inputdir_fx = getInputDirectories(args.dataroot,args.compset,args.experiment)

	# Store correct files (between time boundaries, with correct case name)
	ncfiles = []
	for file in glob.glob(inputdir+"/*.nc"):
		filename = os.path.basename(file)
		if case in filename and ".cam.h0." in filename:
			dt_info = [int(s) for s in filename.split('.')[3].split('-') if s.isdigit()]
			d = datetime(dt_info[0],dt_info[1],dt_info[2],dt_info[3]/3600)
			if d <= dt_bnds[1] and d >= dt_bnds[0]:
				ncfiles.append(file)
	# Exit if no interesting files are found
	if len(ncfiles) == 0:
		print "No files matching request."
		sys.exit(0)
	ncfiles.sort()

	return int(args.fraction), args.experiment, args.compset, args.dataroot, \
		datetime1[:8], datetime2[:8], ncfiles, args.dataroot


def defineOutputDirectory(dataroot,compset,experiment):

	hostname = socket.gethostname()
	case = "bf_%s_%s"%(compset,experiment)
	if hostname == "jollyjumper":
		outputdir = os.path.join(dataroot,'preprocessed',case,'day')
	elif "edison" in hostname:
		outputdir = os.path.join(os.path.dirname(currentpath),'preprocessed',case,'day')

	return outputdir


###--- Main program ---###

if __name__ == "__main__":

	t0 = datetime.now()		# Start chronometer

	##-- Get arguments --##

	fraction, experiment, compset, case, datetime1, datetime2, ncfiles, \
		dataroot = parseArguments()
	case = "bf_%s_%s"%(compset,experiment)

	##-- Variables needed --##

	ref_id = "T"		# reference variable for dimensions and size
	temp_id = "CRM_T"
	pressurf_id = "PS"
	indvarid = "I%d"%fraction
	if fraction < 0:
		indvarid = "all"

	##-- Get input directories --##

	inputdir, inputdir_processed_day, inputdir_processed_1hr, inputdir_results, \
		inputdir_fx = getInputDirectories(dataroot,compset,experiment)

	##-- Output metadata --##

	newvar_id = 'CRM_QVSATENV_%s'%indvarid
	outputdir = defineOutputDirectory(dataroot,compset,experiment)
	outputfile = "%s_day_CESM111-SPCAM20_%s_r1i1p1_%s-%s.nc"%\
		(newvar_id.replace('_','-'),experiment,datetime1,datetime2)	

	##-- Create outputfile --##

	fh_source = Dataset(ncfiles[0],'r')
	lon_source = fh_source.variables['lon']
	lat_source = fh_source.variables['lat']
	lev_source = fh_source.variables['lev']
	ref_dims = fh_source.variables[ref_id].dimensions
	ref_shape = fh_source.variables[ref_id].shape
	temp_shape = fh_source.variables[temp_id].shape
	rootgrp = Dataset(os.path.join(outputdir,outputfile),'w')

		#- Create dimensions and global attributes -#

	rootgrp.createDimension("lon",len(lon_source))
	nlat = len(lat_source)
	lat_slice = slice(nlat/3,2*nlat/3)
	rootgrp.createDimension("lat",nlat/3)
	rootgrp.createDimension("lev",len(fh_source.variables['lev']))
	rootgrp.createDimension("time",None)
	rootgrp.case = str(fh_source.case)
	rootgrp.description = "Daily %s values from CESM111-SPCAM20 in the %s scenario"%\
		(newvar_id,experiment)

		#- Define variables -#

	lon = rootgrp.createVariable("lon","f8",("lon",))
	lon[:] = lon_source[:]
	lat = rootgrp.createVariable("lat","f8",("lat",))
	lat[:] = lat_source[lat_slice]
	lev = rootgrp.createVariable("lev","f8",("lev",))
	lev[:] = fh_source.variables['lev'][:]
	time = rootgrp.createVariable("time","f8",("time",))
	date = rootgrp.createVariable("date","i4",("time",))
	datesec = rootgrp.createVariable("datesec","i4",("time",))
	newvar = rootgrp.createVariable(newvar_id,"f4",ref_dims)

	#- Define variable attributes -#
	lon.long_name = str(lon_source.long_name)
	lat.long_name = str(lat_source.long_name)
	lev.long_name = str(fh_source.variables['lev'].long_name)
	date.long_name = str(fh_source.variables['date'].long_name)
	datesec.long_name = str(fh_source.variables['datesec'].long_name)
	lon.units = str(lon_source.units)
	lat.units = str(lat_source.units)
	lev.units = str(fh_source.variables['lev'].units)
	time.units = str(fh_source.variables['time'].units)
	time.calendar = str(fh_source.variables['time'].calendar)
	if fraction >= 0:
		newvar.long_name = str("Saturation specific humidity from CRM-scale "+\
			"environmental temperature profile and GCM-scale surface pressure "+\
			"in %s rainiest points"%indvarid)
	else:
		newvar.long_name = str("Saturation specific humidity from CRM-scale "+\
			"environmental temperature profile and GCM-scale surface pressure "+\
			"averaged across all CRM points")
	newvar.units = str("kg/kg")

	fh_source.close()


	##-- Read all files --##

		#- Initialize -#

	n = len(ncfiles)
	time_values = np.zeros((n,))
	date_values = np.zeros((n,))
	datesec_values = np.zeros((n,))
	temp_values = np.zeros((n,)+temp_shape[1:-2]+(nlat/3,temp_shape[-1]))
	pressurf_values = np.zeros((n,)+ref_shape[1:-2]+(nlat/3,temp_shape[-1]))

	#- Import required data -#

	for i in range(n):

		ncfile = ncfiles[i]
		
		fh = Dataset(ncfile)		# Open file
		# Copy values
		time_values[i] = fh.variables['time'][:]
		date_values[i] = fh.variables['date'][:]
		datesec_values[i] = fh.variables['datesec'][:]
		temp_values[i] = fh.variables[temp_id][...,lat_slice,:]
		pressurf_values[i] = fh.variables[pressurf_id][...,lat_slice,:]

		fh.close()		# Close file

		#- Import CRM indices -#

	if fraction >= 0:
		ind_values = getVar(indvarid,inputdir_processed_1hr,dates=(datetime1,datetime2))
	else:
		ind_values = np.ones((n,)+temp_shape[2:-2]+(nlat/3,temp_shape[-1]))

	##-- Compute convective-scale qvstar profile --##

		#- Compute the vertical pressure structure at each GCM point -#

	t_pres_profile_0 = datetime.now()
	lev_file = 'lev_fx_CESM111-SPCAM20_allExperiments_r0i0p0.nc'
	computeP = getPressureCoordinateFunctionKeepDimensions(os.path.join(inputdir_fx,lev_file))
	pres_values = computeP(pressurf_values)			# Compute vertical structure
	pres_values = np.moveaxis(pres_values,-1,1)		# Reorder dimensions
	t_pres_profile_1 = datetime.now()
	print "Time to compute p profiles:", t_pres_profile_1-t_pres_profile_0

		#- Extend variables so that they have the same shape -#

	t_reshape_all_0 = datetime.now()
	ntime,nlev,ny,nx,nlat,nlon = temp_values.shape
	# No need to extend temperature values
	# Add crm_y and crm_x dimensions to pressure profiles
	pres_values = np.moveaxis(np.repeat(np.repeat(pres_values[...,np.newaxis],\
												ny,4)[...,np.newaxis],\
												nx,5),\
								[-2,-1],[2,3])
	# Add lev dimension to indices
	ind_values = np.moveaxis(np.repeat(ind_values[...,np.newaxis],nlev,5),-1,1)
	t_reshape_all_1 = datetime.now()
	print "Time to reshape all variables:", t_reshape_all_1-t_reshape_all_0

	# Reverse order or temperature values in the vertical
	temp_values = temp_values[:,::-1,...]

		#- Compute qvstar at each CRM point -#

	t_qvstar_0 = datetime.now()
	newvar_values = saturationSpecificHumidity(temp_values,pres_values)
	t_qvstar_1 = datetime.now()
	print "Time to compute CRM profiles of qvstar:", t_qvstar_1-t_qvstar_0

		#- Compute convective-scale mean qvstar profile -#

	# Replace values at locations of zero ind with NaNs (i.e. outside IXX wet regions)
	newvar_6D_values[ind_values == 0] = np.nan
	# Average in each CRM the values that are not NaNs
	with warnings.catch_warnings():
		warnings.filterwarnings('ignore')
		newvar_4D_values = np.nanmean(np.nanmean(newvar_6D_values,axis=3),axis=2)
	# Now these NaN points correspond to GCM locations where the time slice does
	# not see any of the IXX wet region (checked in the code by looking at the 
	# size of both sets of points)
	# Now average over time the non-NaN values
	newvar_3D_values = np.nanmean(newvar_4D_values,axis=0)
	# Add a time dimension of length 1
	newvar_day_3D_values = newvar_3D_values.reshape((1,)+newvar_3D_values.shape)
	# At this point the NaN values should be locations of zero GCM-scale rainfall
	# (checked on nvciew for the first day of analysis)
	# # Replace the remaining NaN values with zero -- is that necessary?
	# newvar_day_3D_values[np.isnan(omega_values)] = 0.


	##-- Write to output --##

	time[:] = time_values[0]
	date[:] = date_values[0]
	datesec[:] = datesec_values[0]
	newvar[:] = newvar_day_3D_values[:]

	t1 = datetime.now()			# Stop chronometer
	print t1-t0

	rootgrp.close()

	sys.exit(0)


