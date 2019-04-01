################################################################################
#                                                                              #
#    Benjamin Fildier, Nov 2016.                                               #
#                                                                              #
#    Compute the mean profile for 3D variables on the GCM grid for the varids  #
#    given in input, and save them in a data frame.                            #
#                                                                              #
################################################################################


###--- Modules ---###
import sys,os
from netCDF4 import Dataset
import numpy as np
import socket
from argparse import ArgumentParser
from pandas import DataFrame
import csv

## Own functions
currentpath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,currentpath[:currentpath.rfind('/')+1]+'functions')

from importingData import *
from statistics import *


###--- Functions ---###

def parseArguments():

	parser = ArgumentParser(description="Compute mean profiles for a sequence of 3D variables")
	parser.add_argument('varids',nargs="+",
		help="Sequence of 3D variables defined on GCM grid to analyze")
	parser.add_argument('-e','--experiment',default="piControl",
		choices=['piControl','abrupt4xCO2'],
		help="Physical experiment or scenario.")
	parser.add_argument('-s','--subset',default="tropics",
		choices=['tropics','ocean','land','mfzero'],
		help="Reduce the domain to a subset of points.")
	parser.add_argument('-d','--directory',required=True,
		help="Input directory root")
	parser.add_argument('-c','--compset',choices=['FSPCAMm_AMIP','FAMIPC5'],
		default="FSPCAMm_AMIP",help="CESM compset used.")
	args = parser.parse_args()

	return args.varids, args.directory, args.compset, args.experiment, args.subset

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

def getSubset(subsetName,inputdir_fx,experiment=None):

	if subsetName == "ocean":
		# Define subset of points over ocean
		landfile = 'landmask_fx_CESM111-SPCAM20_allExperiments_r0i0p0.nc'
		fh_landmask = Dataset(os.path.join(inputdir_fx,landfile))
		landmask = fh_landmask.variables['landmask']
		subset_pts = np.logical_not(landmask)
		fh_landmask.close()
	elif subsetName == "land":
		# Define subset of points over land
		landfile = 'landmask_fx_CESM111-SPCAM20_allExperiments_r0i0p0.nc'
		fh_landmask = Dataset(os.path.join(inputdir_fx,landfile))
		landmask = fh_landmask.variables['landmask']
		subset_pts = np.array(landmask[:],dtype=bool)
		fh_landmask.close()
	elif subsetName == "tropics":
		subset_pts = None
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

	varids, dataroot, compset, experiment, subset = parseArguments()

	##-- Set environment --##

	inputdir, inputdir_fx = getInputDirectories(dataroot,compset,experiment)
	outputdir = os.path.join(os.path.dirname(currentpath),'results')
	outputfile = "var3D_meanProfile_%s_%s_%s.csv"%(compset,experiment,subset)

	##-- Define subset if necessary --##

	subset_pts = getSubset(subset,inputdir_fx,experiment=experiment)

	##-- Compute means --##

	var_meanProfile = {}
	removeNanProfiles = True

	for varid in varids:

		print "  Target variable:", varid

		##-- Import data --##

		var = getVar(varid,inputdir)
		if var is None:
			print varid, "ignored."
			continue
		elif len(var.shape) == 3:
			print varid, "is 2D, ignored."
			continue

		##-- Compute means --##

		meanProfile = computeMean(var,subset=subset_pts)

		## Save if vertical profile is not full of nans, otherwise remove it
		if np.isnan(meanProfile).sum() < meanProfile.size or not removeNanProfiles:
			var_meanProfile[varid] = meanProfile
		else:
			print "    Mean profile of %s is Nan, ignored."%varid

	##-- Turn results into dataframes and save to CSV files --##

	## Convert to dataframe
	df_mean = DataFrame(var_meanProfile)

	## Save to CSV
	df_mean.to_csv(os.path.join(outputdir,outputfile))





