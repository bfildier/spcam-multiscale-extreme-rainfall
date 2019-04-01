################################################################################
#                                                                              #
#    Benjamin Fildier, March 2017.                                             #
#                                                                              #
#    Compute GCM indices corresponding to a given percentile of the rainfall   #
#    distribution, for a given reference PR_ID, a given percentile value Q_ID  #
#    and the experiment, compset and subset considered (tropics, land, ocean). #
#                                                                              #
#    Arguments   : PR_ID Q_ID_bnds inputdir subset experiment compset datetimes#
#                                                                              #
################################################################################

###--- Modules ---###
import sys,os
from netCDF4 import Dataset
import numpy as np
import socket
from argparse import ArgumentParser
import csv
import pandas as pd

## Own functions
currentpath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,currentpath[:currentpath.rfind('/')+1]+'functions')

from importingData import *
from statistics import *


###--- Functions ---###

def parseArguments():

	parser = ArgumentParser(description="Extract GCM points corresponding to a given quantile of precipitation.")
	parser.add_argument('-p','--pr_id',required=True,
		help='Reference precipitation ID.')
	parser.add_argument('-q','--qRanks',nargs="+",required=True,help='quantile')
	parser.add_argument('-d','--directory',required=True,
		help="Input directory root.")
	parser.add_argument('-e','--experiment',default="piControl",
		choices=['piControl','abrupt4xCO2'],
		help="Physical experiment or scenario.")
	parser.add_argument('-s','--subset',default="tropics",
		choices=['tropics','ocean','land','mfzero'],
		help="Reduce the domain to a subset of points.")
	parser.add_argument('-c','--compset',choices=['FSPCAMm_AMIP','FAMIPC5'],
		default="FSPCAMm_AMIP",help="CESM compset used.")
	parser.add_argument('-dt','--datetimes',nargs=2,required=True,
		help="Boundary datetime values in YYYY-MM-DD or YYYY-MM-DD-SSSSS format.")
	args = parser.parse_args()

	# Get time boundaries
	dt_bnds = []
	for tstamp in args.datetimes:
		dt_info = [int(s) for s in tstamp.split('-') if s.isdigit()]
		dt_bnds.append(datetime(dt_info[0],dt_info[1],dt_info[2]))
	# Change date format to CMIP5 standard
	datetime1 = dt_bnds[0].isoformat().replace('-','').replace('T','').replace(':','')[:-2]
	datetime2 = dt_bnds[1].isoformat().replace('-','').replace('T','').replace(':','')[:-2]

	return args.pr_id, args.qRanks, args.directory, args.compset, args.experiment, \
		args.subset, (datetime1,datetime2)

def getInputDirectories(dataroot,compset,experiment):

	hostname = socket.gethostname()
	case = "bf_%s_%s"%(compset,experiment)
	if hostname == "jollyjumper":
		inputdir = os.path.join(dataroot,'preprocessed',case,'day')
	elif "edison" in hostname or "cori" in hostname:	
		inputdir = os.path.join(os.path.dirname(currentpath),'preprocessed',
			case,'day')
	inputdir_fx = os.path.join(os.path.dirname(os.path.dirname(inputdir)),
		'allExperiments/fx')

	return inputdir, inputdir_fx

def getPercentileFile(referencepath,compset,experiment,subset):

	precdir = os.path.join(os.path.dirname(referencepath),'results')
	precfile = "pr_quantile_IL_%s_%s_%s.csv"%(compset,experiment,subset)

	return precdir,precfile

def getSubset(subsetName,inputdir_fx,experiment=None):

	landfile = 'landmask_fx_CESM111-SPCAM20_allExperiments_r0i0p0.nc'
	fh_landmask = Dataset(os.path.join(inputdir_fx,landfile))
	landmask = fh_landmask.variables['landmask']

	if subsetName == "ocean":
		# Define subset of points over ocean
		subset_pts = np.logical_not(landmask)
		fh_landmask.close()
	elif subsetName == "land":
		# Define subset of points over land
		subset_pts = np.array(landmask[:],dtype=bool)
		fh_landmask.close()
	elif subsetName == "tropics":
		subset_pts = np.ones(landmask.shape,dtype=bool)
	elif subsetName == "mfzero":
		# mfzero_file = 'MFZERO_day_CESM111-SPCAM20_'+experiment+'_r1i1p1_18500501-18500502.nc'
		mfzero_file = 'MFZERO_day_CESM111-SPCAM20_'+experiment+'_r1i1p1_18500501-18510430.nc'
		fh_mfzero = Dataset(os.path.join(os.path.dirname(inputdir_fx),'day',mfzero_file))
		mfzero = fh_mfzero.variables['MFZERO']
		subset_pts = np.array(mfzero[:],dtype=bool)
		fh_mfzero.close()

	return subset_pts


###---- Main program ---###

if __name__ == "__main__":

	##-- Get arguments --##

	pr_id, qRanks, dataroot, compset, experiment, subsetName, datetimes =  parseArguments()
	if compset == 'FAMIPC5':
		model = 'CESM111-CAM5'
	else:
		model = 'CESM111-SPCAM20'

	##-- Set environment --##

	inputdir, inputdir_fx = getInputDirectories(dataroot,compset,experiment)
	percdir, percfile = getPercentileFile(currentpath,compset,experiment,subsetName)
	outputdir = inputdir

	##-- Get bounds of precipitation range

	df_perc = pd.read_csv(os.path.join(percdir,percfile))
	Qs = np.array([float(Q) for Q in df_perc['Q_IL']])
	Qmin = float(qRanks[0])
	i = np.argmin(np.absolute(Qs-Qmin))
	if len(qRanks) > 1:
		Qmax = float(qRanks[1])
		j = np.argmin(np.absolute(Qs-Qmax))
	elif i < len(Qs)-1:
		j = i+1
	else:	
		j = 0
	pr_min = df_perc[pr_id][i]
	pr_max = df_perc[pr_id][j]

	##-- Get daily precipitation values

	pr_vals = getVar(pr_id,inputdir,
		dates=(datetimes[0].replace('-',''),datetimes[1].replace('-','')))*86400000 # convert to m/day

	##-- Define subset if necessary --##

	subset = getSubset(subsetName,inputdir_fx,experiment=experiment)

	##-- Define new var id and outputfile --##

	newvarid = "J%s_%s"%(string.join(("%2.4f"%Qs[i]).split('.'),''),pr_id.replace('-','_'))
	outputfile = "%s_day-%s_%s_%s_r1i1p1_%s-%s.nc"%(newvarid.replace('_','-'),subsetName,\
	model,experiment,datetimes[0][:8],datetimes[1][:8])

	##-- Create outputfile --##

	#- Get original dimensions -#
	ncfiles = getInputfiles(pr_id,inputdir,
		dates=(datetimes[0].replace('-',''),datetimes[1].replace('-','')))
	fh_source = Dataset(ncfiles[0],'r')
	lon_source = fh_source.variables['lon']
	lat_source = fh_source.variables['lat']
	vardims = fh_source.variables[pr_id].dimensions
	varshape = fh_source.variables[pr_id].shape

	#- Create file -#
	rootgrp = Dataset(os.path.join(outputdir,outputfile),'w')
	#- Create dimensions and global attributes -#
	nlon = len(lon_source)
	rootgrp.createDimension("lon",nlon)
	nlat = len(lat_source)
	rootgrp.createDimension("lat",nlat)
	rootgrp.createDimension("time",None)
	rootgrp.case = str(fh_source.case)
	rootgrp.description = ("GCM points with %s values at their %sth "+\
	"percentile for the %s subdomain and %s experiment.")%\
	(pr_id,qRanks[0],subsetName,experiment)
	
	#- Define variables -#
	lon = rootgrp.createVariable("lon","f8",("lon",))
	lon[:] = lon_source[:]
	lat = rootgrp.createVariable("lat","f8",("lat",))
	lat[:] = lat_source[:]
	time = rootgrp.createVariable("time","f8",("time",))
	time[:] = fh_source.variables['time'][:]
	date = rootgrp.createVariable("date","i4",("time",))
	date[:] = fh_source.variables['date'][:]
	datesec = rootgrp.createVariable("datesec","i4",("time",))
	datesec[:] = fh_source.variables['datesec'][:]
	newvar = rootgrp.createVariable(newvarid,"u1",("time","lat","lon"))

	#- Define variable attributes -#
	lon.long_name = str(lon_source.long_name)
	lat.long_name = str(lat_source.long_name)
	date.long_name = str(fh_source.variables['date'].long_name)
	datesec.long_name = str(fh_source.variables['datesec'].long_name)
	lon.units = str(lon_source.units)
	lat.units = str(lat_source.units)
	time.units = str(fh_source.variables['time'].units)
	time.calendar = str(fh_source.variables['time'].calendar)
	newvar.long_name = "Indices of %s %s percentile"%(pr_id,qRanks[0])
	newvar.units = ""

	fh_source.close()

	##-- Compute mask of correct indices, intersect with subset and store to output --##

	# Location of percentiles
	if pr_max > pr_min:
		loc_perc = np.array(np.logical_and(pr_vals >= pr_min, pr_vals < pr_max),dtype=bool)
	else:
		loc_perc = np.array(pr_vals >= pr_min,dtype=bool)
	# Intersect with subset
	loc_perc_subset = np.logical_and(loc_perc,subset)
	loc_perc_subset = np.array(loc_perc_subset,dtype=int)
	# Store values
	newvar[:] = loc_perc_subset[:]

	rootgrp.close()

	sys.exit(0)















