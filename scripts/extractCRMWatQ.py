################################################################################
#                                                                              #
#    Benjamin Fildier, July 2017.                                              #
#                                                                              #
#    Extract CRM-W-IXX profiles at GCM locations of CRM_PREC_IXX quantile $Q_ID#
#    within the fraction area corresponding to the IXX valid indices. For a    #
#    given experiment, compset and subset considered (tropics, land, ocean).   #
#                                                                              #
#    Arguments   : PR_ID Q_ID dataroot subset experiment compset date          #
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
import datetime as dt
from datetime import datetime,timedelta

## Own functions
currentpath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,currentpath[:currentpath.rfind('/')+1]+'functions')

from importingData import *
from statistics import *
from thermodynamics import airDensity,g
from extremeScaling import verticalPressureIntegral,verticalPressureIntegralProduct


###--- Functions ---###

def parseArguments():

	parser = ArgumentParser(description="Extract GCM points corresponding to a given quantile of precipitation.")
	parser.add_argument('-p','--pr_id',required=True,
		help='Reference precipitation ID.')
	parser.add_argument('-a','--area_id',required=True,
		help='Fraction area defining ID')
	parser.add_argument('-q','--q_id',required=True,help='quantile')
	parser.add_argument('-d','--dataroot',required=True,
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

	return args.pr_id.replace('-','_'), args.area_id, args.q_id, args.dataroot, args.compset, args.experiment, \
		args.subset, datetime1[:8], datetime2[:8]

def getInputDirectories(dataroot,compset,experiment):

	hostname = socket.gethostname()
	case = "bf_%s_%s"%(compset,experiment)
	if hostname == "jollyjumper":
		inputdir = os.path.join(dataroot,"simulations",case)
		inputdir_processed_day = os.path.join(dataroot,'preprocessed',case,'day')
		inputdir_processed_1hr = os.path.join(dataroot,'preprocessed',case,'1hr')
	elif "edison" in hostname or "cori" in hostname:
		inputdir = os.path.join(dataroot,'archive',case,"atm/hist")
		inputdir_processed_day = os.path.join(os.path.dirname(currentpath),'preprocessed',
			case,'day')
		inputdir_processed_1hr = os.path.join(os.path.dirname(currentpath),'preprocessed',
			case,'1hr')
	inputdir_results = os.path.join(os.path.dirname(currentpath),'results')
	inputdir_fx = os.path.join(dataroot,'preprocessed/allExperiments/fx')

	return inputdir, inputdir_processed_day, inputdir_processed_1hr,\
		inputdir_results, inputdir_fx

def defineOutputDirectory(dataroot,compset,experiment):

	hostname = socket.gethostname()
	case = "bf_%s_%s"%(compset,experiment)
	if hostname == "jollyjumper":
		outputdir = os.path.join(os.path.dirname(os.path.dirname(inputdir)),
			'preprocessed',case,'day')
	elif "edison" in hostname or "cori" in hostname:
		outputdir = os.path.join(os.path.dirname(currentpath),'preprocessed',case,'day')

	return outputdir

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


###--- Main program ---###

if __name__ == "__main__":

	t0 = datetime.now()		# Start chronometer

##-- Get arguments --##

	pr_id, area_id, q_id, dataroot, compset, experiment, subsetName, datetime1, datetime2 =  parseArguments()
	if compset == 'FAMIPC5':
		model = 'CESM111-CAM5'
	else:
		model = 'CESM111-SPCAM20'

##-- Set directories --##

	inputdir, inputdir_processed_day, inputdir_processed_1hr, inputdir_results,\
		inputdir_fx = getInputDirectories(dataroot,compset,experiment)
	outputdir = defineOutputDirectory(dataroot,compset,experiment)

##-- Load required physical variables -- ##

	w_id = "CRM_W_%s"%area_id
	# Load all variables
	w_vals = getVar(w_id,inputdir_processed_day,dates=(datetime1,datetime2))

##-- Load required index information --##

	# Index J of quantile
	percfile = "pr_quantile_IL_%s_%s_%s.csv"%(compset,experiment,subsetName)
	df_perc = pd.read_csv(os.path.join(inputdir_results,percfile))
	Qs = np.array([float(Q) for Q in df_perc['Q_IL']])
	J = np.argmin(np.absolute(Qs-float(q_id)))
	Q_id = ("%2.4f"%Qs[J]).replace('.','')
	J_id = "J_%s"%(pr_id)

	# J-prid values
	J_vals = getVar(J_id,inputdir_processed_day,dates=(datetime1,datetime2))[0]
	J_valid = np.array(J_vals == J,dtype=bool)

##-- Load subset --##

	subset = getSubset(subsetName,inputdir_fx,experiment=experiment)

##-- Define variable, new file name --##

	newvarid = "CRM_W_%s_%s_%s"%(area_id,pr_id.replace('-','_'),Q_id)
	outputfile = '%s_day-%s_%s_%s_r1i1p1_%s_%s.nc'%(newvarid.replace('_','-'),subsetName,model,\
		experiment,datetime1,datetime2)

##-- Create and initialize new file --##

	lev_file = 'lev_fx_CESM111-SPCAM20_allExperiments_r0i0p0.nc'
	fh_source = Dataset(os.path.join(inputdir_fx,lev_file),'r')
	lev_source = fh_source.variables['lev']
	nlev = len(lev_source)

	# Compute dimension Nranks	
	Nranks = J_valid.sum()
	
	# Create outputfile and dimensions
	rootgrp = Dataset(os.path.join(outputdir,outputfile),'w')
	rootgrp.createDimension("rank",Nranks)
	rootgrp.createDimension("lev",nlev)
	rootgrp.createDimension("time",None)
	rootgrp.case = str(fh_source.case)
	rootgrp.description = ("Convective-scale verical velocity profiles for %s wet area "+\
		"at location of %s's %2.4fth percentiles.")%(area_id,pr_id,Qs[J])

	# Define variables
	crm_w = rootgrp.createVariable(newvarid,"f8",("time","lev","rank"))
	date = rootgrp.createVariable("date","i4",("time",))
	crm_w.long_name = "Convective-scale vertical velocities"
	date.long_name = "current date (YYYYMMDD)"
	crm_w.units = 'm/s'
	date[:] = datetime1

##-- Get convective-scale W profile at relevant points --##
	
		#- Extend variables so that they have the same shape -#

	# Add lev dimension to J_valid
	J_valid = np.moveaxis(np.repeat(J_valid[...,np.newaxis],nlev,len(J_valid.shape)),[-1],[0])
	# Add time dimension to J_valid
	J_valid = J_valid.reshape((1,)+J_valid.shape)

		#- Extract profiles at quantile -#

	w_vals_at_Q = w_vals[J_valid]
	w_vals_at_Q = np.reshape(w_vals_at_Q,(1,nlev,Nranks))

	crm_w[:] = w_vals_at_Q

	t1 = datetime.now()			# Stop chronometer
	print t1-t0
	
	rootgrp.close()

	sys.exit(0)
