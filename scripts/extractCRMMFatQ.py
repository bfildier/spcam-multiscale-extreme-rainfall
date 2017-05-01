################################################################################
#                                                                              #
#    Benjamin Fildier, March 2017.                                             #
#                                                                              #
#    Compute and extract CRM-level column integrated mass flux values at GCM   #
#    locations of CRM_PREC_IXX quantile $Q_ID within the fraction area         #
#    corresponding to the IXX valid indices. For a given experiment, compset   #
#    and subset considered (tropics, land, ocean).                             #
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
		args.subset, datetime1, datetime2

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

def defineOutputDirectory(dataroot,compset,experiment):

	hostname = socket.gethostname()
	case = "bf_%s_%s"%(compset,experiment)
	if hostname == "jollyjumper":
		outputdir = os.path.join(os.path.dirname(os.path.dirname(inputdir)),
			'preprocessed',case,'1hr')
	elif "edison" in hostname:
		outputdir = os.path.join(os.path.dirname(currentpath),'preprocessed',case,'1hr')

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
	
	w_id = "CRM_W"
	spechum_id = 'Q'
	temp_id = 'T'
	pressurf_id = 'PS'
	pw_id = 'CRM_QPC'
	# Load all variables
	w_vals = getSimulationValues(w_id,inputdir,dates=(datetime1,datetime2))
	spechum_vals = getSimulationValues(spechum_id,inputdir,dates=(datetime1,datetime2))
	temp_vals = getSimulationValues(temp_id,inputdir,dates=(datetime1,datetime2))
	pressurf_vals = getSimulationValues(pressurf_id,inputdir,dates=(datetime1,datetime2))
	pw_vals = getSimulationValues(pw_id,inputdir,dates=(datetime1,datetime2))
	# Reduce latitude dimension
	nlat = pressurf_vals.shape[-2]
	lat_slice = slice(nlat/3,2*nlat/3)
	w_vals = w_vals[:,:,:,:,lat_slice,:]
	spechum_vals = spechum_vals[:,:,lat_slice,:]
	temp_vals = temp_vals[:,:,lat_slice,:]
	pressurf_vals = pressurf_vals[:,lat_slice,:]
	pw_vals = pw_vals[:,:,:,:,lat_slice,:]

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

	# IXX indices
	I_vals = getVar(area_id,inputdir_processed_1hr,dates=(datetime1,datetime2))

##-- Load computeP and subset --##

	lev_file = 'lev_fx_CESM111-SPCAM20_allExperiments_r0i0p0.nc'
	computeP = getPressureCoordinateFunction(os.path.join(inputdir_fx,lev_file))
	subset = getSubset(subsetName,inputdir_fx,experiment=experiment)

##-- Define variable, new file name --##

	newvarid = "CRM_MF_%s_%s_%s"%(area_id,pr_id.replace('-','_'),Q_id)
	outputfile = '%s_1hr-%s_%s_%s_r1i1p1_%s_%s.nc'%(newvarid.replace('_','-'),subsetName,model,\
		experiment,datetime1,datetime2)

##-- Create and initialize new file --##

	# Compute dimension Nranks	
	I_vals_summed = I_vals.sum(axis=(0,1,2)) # Sum along time, crm_ny and crm_nx dimensions
	Nranks = np.diagonal(np.dot(I_vals_summed,J_valid.transpose())).sum()
	
	# Create outputfile and dimensions
	rootgrp = Dataset(os.path.join(outputdir,outputfile),'w')
	rootgrp.createDimension("time",None)
	rootgrp.createDimension("rank",Nranks)
	rootgrp.description = "Mass flux values at CRM points contributing to %s's %2.4fth percentile"%\
		(pr_id,Qs[J])

	# Define variables
	is_new_cell = rootgrp.createVariable("isNewCell","u1",("time","rank"))
	crm_mf = rootgrp.createVariable(newvarid,"f8",("time","rank"))
	crm_mfup = rootgrp.createVariable(newvarid.replace('MF','MFUP'),"f8",("time","rank"))
	crm_mfdn = rootgrp.createVariable(newvarid.replace('MF','MFDN'),"f8",("time","rank"))
	crm_mf_pw = rootgrp.createVariable(newvarid.replace('MF','MF_PW'),"f8",("time","rank"))
	crm_mfup_pw = rootgrp.createVariable(newvarid.replace('MF','MFUP_PW'),"f8",("time","rank"))
	crm_mfdn_pw = rootgrp.createVariable(newvarid.replace('MF','MFDN_PW'),"f8",("time","rank"))
	crm_pw = rootgrp.createVariable(newvarid.replace('MF','PW'),"f8",("time","rank"))
	date = rootgrp.createVariable("date","i4",("time",))
	is_new_cell.long_name = "Start of new GCM cell"
	crm_mf.long_name = "CRM mass flux"
	crm_mfup.long_name = "CRM mass flux up"
	crm_mfdn.long_name = "CRM mass flux down"
	crm_mf_pw.long_name = "Precipitable water-weighted CRM mass flux"
	crm_mfup_pw.long_name = "Precipitable water-weighted CRM mass flux up"
	crm_mfdn_pw.long_name = "Precipitable water-weighted CRM mass flux down"
	crm_pw.long_name = "CRM precipitable water"
	date.long_name = "current date (YYYYMMDD)"
	is_new_cell.units = ''
	crm_mf.units = 'Pa/s'
	date[:] = datetime1[:8]

##-- Compute CRM mass flux at relevant points --##

	# Iinitialize rank
	k = 0
	is_new_cell[0,:] = np.zeros((1,Nranks))
	# Loop over GCM points that match the percentile of interest
	for ilat,ilon in zip(*np.where(J_valid)):
		# Mark that we entered a new GCM cell
		is_new_cell[0,k] = 1
		# Get the CRM points that correspond to the convective event
		I_valid = np.array(I_vals[...,ilat,ilon] == 1,dtype=bool).squeeze()
		for itime,ix in zip(*np.where(I_valid)):
			# Compute and store CRM mass flux at current CRM point
			pressurf = pressurf_vals[itime,ilat,ilon]
			pres = computeP(pressurf)
			spechum = spechum_vals[itime,:,ilat,ilon]
			temp = temp_vals[itime,:,ilat,ilon]
			rho = airDensity(temp,pres,spechum)
			wspeed = w_vals[itime,...,ix,ilat,ilon].squeeze()
			pw = pw_vals[itime,...,ix,ilat,ilon].squeeze()
			omega = np.multiply(wspeed,-g*rho)
			omega_up = omega.copy()
			omega_dn = omega.copy()
			omega_up[omega >= 0] = 0
			omega_dn[omega <= 0] = 0
			crm_mf[0,k] = verticalPressureIntegral(pres,omega) / verticalPressureIntegral(pres)
			crm_mfup[0,k] = verticalPressureIntegral(pres,omega_up) / verticalPressureIntegral(pres)
			crm_mfdn[0,k] = verticalPressureIntegral(pres,omega_dn) / verticalPressureIntegral(pres)
			crm_mf_pw[0,k] = verticalPressureIntegralProduct(pres,omega,pw) / verticalPressureIntegral(pres,pw)
			crm_mfup_pw[0,k] = verticalPressureIntegralProduct(pres,omega_up,pw) / verticalPressureIntegral(pres,pw)
			crm_mfdn_pw[0,k] = verticalPressureIntegralProduct(pres,omega_dn,pw) / verticalPressureIntegral(pres,pw)
			crm_pw[0,k] = verticalPressureIntegral(pres,pw) / verticalPressureIntegral(pres)
			# increment current rank
			k += 1
	
	rootgrp.close()

	sys.exit(0)
