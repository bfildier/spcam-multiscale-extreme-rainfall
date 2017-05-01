################################################################################
#                                                                              #
#    Benjamin Fildier, Nov 2016.                                               #
#                                                                              #
#    Compute the precipitation quantiles on an inverse-log axis of quantile    #
#    ranks for the set of precipitation variables given in input, and save     #
#    them in a data frame.                                                     #
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

	parser = ArgumentParser(description="Compute high quantiles and rank for a sequence of precipitation variables")
	parser.add_argument('pr_ids',nargs="+",
		help="Sequence of precipitation variables IDs to analyze")
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

	return args.pr_ids, args.directory, args.compset, args.experiment, args.subset

def getInputDirectories(dataroot,compset,experiment):

	hostname = socket.gethostname()
	case = "bf_%s_%s"%(compset,experiment)
	if hostname == "jollyjumper":
		inputdir = os.path.join(dataroot,'preprocessed',case,'day')
	elif "edison" in hostname:	
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

	pr_ids, dataroot, compset, experiment, subset = parseArguments()

	##-- Set environment --##

	inputdir, inputdir_fx = getInputDirectories(dataroot,compset,experiment)
	outputdir = os.path.join(os.path.dirname(currentpath),'results')
	outputfile_quantiles_IL = "pr_quantile_IL_%s_%s_%s.csv"%(compset,experiment,subset)
	outputfile_density_IL = "pr_density_IL_%s_%s_%s.csv"%(compset,experiment,subset)
	# outputfile_mids_L = "pr_mids_L_%s_%s_%s.csv"%(compset,experiment,subset)
	# outputfile_density_L = "pr_density_L_%s_%s_%s.csv"%(compset,experiment,subset)

	##-- Define subset if necessary --##

	subset_pts = getSubset(subset,inputdir_fx,experiment=experiment)

	##-- Compute precipitation quantiles --##

	pr_quantiles_IL = {}
	pr_density_IL = {}
	# pr_mids_L = {}
	# pr_density_L = {}

	for pr_id in pr_ids:

		##-- Import data --##

		pr = getVar(pr_id,inputdir)
		if pr is not None:
			pr *= 86400*1000    # Convert to mm/day
		else:
			continue

		##-- Compute statistics --##

		maxQ, Q_IL, pr_c_IL, pr_e_IL, pr_c_L, pr_e_L, pr_Hd_L, pr_Hd_IL =\
			getStatistics(pr,subset=subset_pts,n_unit_min=3)

		##-- Store results in dictionaries

		pr_quantiles_IL[pr_id] = pr_c_IL
		pr_density_IL[pr_id] = pr_Hd_IL
		# pr_mids_L[pr_id] = pr_c_L
		# pr_density_L[pr_id] = pr_Hd_L

	##-- Turn results into dataframes and save to CSV files --##

	## Add Q coordinate
	pr_quantiles_IL['Q_IL'] = np.array([str(Q.round(4)) for Q in Q_IL])
	pr_density_IL['Q_IL'] = np.array([str(Q.round(4)) for Q in Q_IL])

	## Convert to dataframes
	df_quantiles_IL = DataFrame(pr_quantiles_IL)
	df_density_IL = DataFrame(pr_density_IL)
	# df_mids_L = DataFrame(pr_mids_L)
	# df_density_L = DataFrame(pr_density_L)

	## Save to csv
	df_quantiles_IL.to_csv(os.path.join(outputdir,outputfile_quantiles_IL))
	df_density_IL.to_csv(os.path.join(outputdir,outputfile_density_IL))
	# df_mids_L.to_csv(os.path.join(outputdir,outputfile_mids_L))
	# df_density_L.to_csv(os.path.join(outputdir,outputfile_density_L))





