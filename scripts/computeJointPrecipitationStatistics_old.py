################################################################################
#                                                                              #
#    Benjamin Fildier, Nov 2016.                                               #
#                                                                              #
#    Compute joint precipitation statistics on inverse-log axes of quantile    #
#    ranks for paris of precipitation variables given in input, and save       #
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
	parser.add_argument('-p1','--pr_ids_1',nargs="+",required=True,
		help="First sequence of precipitation variables IDs")
	parser.add_argument('-p2','--pr_ids_2',nargs="+",required=True,
		help="Second sequence of precipitation variables IDs")
	parser.add_argument('-e','--experiment',default="piControl",
		choices=['piControl','abrupt4xCO2'],
		help="Physical experiment or scenario.")
	parser.add_argument('-s','--subset',default="tropics",
		choices=['tropics','ocean','land'],
		help="Reduce the domain to a subset of points.")
	parser.add_argument('-d','--directory',required=True,
		help="Input directory root")
	parser.add_argument('-c','--compset',choices=['FSPCAMm_AMIP','FAMIPC5'],
		default="FSPCAMm_AMIP",help="CESM compset used.")
	args = parser.parse_args()

	return args.pr_ids_1, args.pr_ids_2, args.directory, args.compset, args.experiment, args.subset

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

def getSubset(subsetName,inputdir_fx):

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

	return subset_pts


###--- Main program ---###

if __name__ == "__main__":

	##-- Get arguments --##

	pr_ids_1, pr_ids_2, dataroot, compset, experiment, subset = parseArguments()

	##-- Set environment --##

	inputdir, inputdir_fx = getInputDirectories(dataroot,compset,experiment)
	outputdir = os.path.join(os.path.dirname(currentpath),'results')
	outputfile_density = "pr_jointDensity_IL_%s_%s_%s.csv"%(compset,experiment,subset)

	##-- Define subset if necessary --##

	subset_pts = getSubset(subset,inputdir_fx)

	##-- Compute precipitation quantiles --##

	pr_jointDensity_IL = {}

	for pr_id_1 in pr_ids_1:

		pr_jointDensity_IL[pr_id_1] = {}

		for pr_id_2 in pr_ids_2:

			##-- Import data --##

			pr1 = getVar(pr_id_1,inputdir)
			pr2 = getVar(pr_id_2,inputdir)
			if pr1 is not None:
				pr1 *= 86400*1000    # Convert to mm/day
				if pr2 is not None:
					pr2 *= 86400*1000    # Convert to mm/day

					##-- Compute statistics --##

					jointstats = getJointStatistics(pr1,pr2,subset=subset_pts,n_unit_min=3)
					if jointstats is None:
						continue
					else:
						maxQ1, maxQ2, Q1_IL, Q2_IL, pr1_e_IL, pr1_c_IL, pr2_e_IL, pr2_c_IL, pr_H2Dd_PI = jointstats

					##-- Store results in dictionaries

					pr_jointDensity_IL[pr_id_1][pr_id_2] = np.array_str(np.array(pr_H2Dd_PI))

	##-- Turn results into dataframes and save to CSV files --##

	## Add Q coordinate
	# pr_jointDensity_IL['Q_IL'] = Q_IL

	## Convert to dataframes
	df_density_IL = DataFrame(pr_jointDensity_IL)
	# df_density_IL = df_density_IL.set_index(df_density_IL.columns[0])

	## Save to csv
	df_density_IL.to_csv(os.path.join(outputdir,outputfile_density))










