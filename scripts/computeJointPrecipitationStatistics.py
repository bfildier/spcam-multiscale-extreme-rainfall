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
from pandas import DataFrame,MultiIndex
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
		choices=['tropics','ocean','land','mfzero'],
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

def stackInPridDimensions(var):

	# print np.array(var).shape
	
	if np.array(var).shape[0] == 1:
		var_stacked_once = np.array(var[0])
	else:
		var_stacked_once = np.vstack(var)

	# print var_stacked_once.shape

	if np.array(var_stacked_once).shape[0] == 1:
		var_stacked = var_stacked_once[0][0]
	else:
		var_stacked = np.hstack(var_stacked_once)[0]

	# print var_stacked.shape

	return var_stacked

###--- Main program ---###

if __name__ == "__main__":

	##-- Get arguments --##

	pr_ids_1, pr_ids_2, dataroot, compset, experiment, subset = parseArguments()

	##-- Set environment --##

	inputdir, inputdir_fx = getInputDirectories(dataroot,compset,experiment)
	outputdir = os.path.join(os.path.dirname(currentpath),'results')
	outputfile_density = "pr_jointDensity_IL_%s_%s_%s.csv"%(compset,experiment,subset)

	##-- Define subset if necessary --##

	subset_pts = getSubset(subset,inputdir_fx,experiment=experiment)

	##-- Compute precipitation quantiles --##

	pr_jointDensity_IL = []
	pr1_to_remove = []
	pr2_to_remove = []

	for pr_id_1 in pr_ids_1:
		# print pr_id_1,

		pr1 = getVar(pr_id_1,inputdir)
		if pr1 is not None:
			pr1 *= 86400*1000    # Convert to mm/day

			pr_jointDensity_IL.append([])
			i = np.argmax(np.array(pr_ids_1) == pr_id_1) - len(pr1_to_remove)

			for pr_id_2 in pr_ids_2:
				# print pr_id_2,

				##-- Import data --##
				
				pr2 = getVar(pr_id_2,inputdir)
				
				if pr2 is not None:
					pr2 *= 86400*1000    # Convert to mm/day

					pr_jointDensity_IL[i].append([])
					j = np.argmax(np.array(pr_ids_2) == pr_id_2) - len(pr2_to_remove)

					##-- Compute statistics --##

					jointstats = getJointStatistics(pr1,pr2,subset=subset_pts,n_unit_min=3)
					if jointstats is None:
						continue
					else:
						maxQ1, maxQ2, Q1_IL, Q2_IL, pr1_e_IL, pr1_c_IL, pr2_e_IL, pr2_c_IL, pr_H2Dd_PI = jointstats

					##-- Store results in list

					H = pr_H2Dd_PI
					# print H.shape
					# print np.array(H.tolist()).shape
					pr_jointDensity_IL[i][j].append(H.tolist())

				else:
					pr2_to_remove.append(pr_id_2)
		else:
			pr1_to_remove.append(pr_id_1)

	## Remove incorrect variables
	for pr_id_1 in pr1_to_remove:
		pr_ids_1.remove(pr_id_1)
	for pr_id_2 in pr2_to_remove:
		pr_ids_2.remove(pr_id_2)

	print "Remaining precipitation variables x-axis:", pr_ids_1
	print "Remaining precipitation variables y-axis:", pr_ids_2

	##-- Turn results into dataframes and save to CSV files --##

	## Define Q values for the first variable
	Q_ids = [str(Q.round(4)) for Q in Q1_IL]

	## Stack twice to get a 1D array of lists
	var_stacked = stackInPridDimensions(pr_jointDensity_IL)

	## Define multidimensional indices for dataframe
	prs1 = [[pr_id_1]*len(pr_ids_2)*len(Q_ids) for pr_id_1 in pr_ids_1]
	prs1 = np.hstack(prs1)
	prs2 = [[pr_id_2]*len(Q_ids) for pr_id_2 in pr_ids_2]
	prs2 = np.hstack(prs2)
	prs2 = np.hstack([prs2]*len(pr_ids_1))
	Qs = Q_ids*len(pr_ids_1)*len(pr_ids_2)

	ind_tuples = zip(prs1,prs2,Qs)
	m_index = MultiIndex.from_tuples(ind_tuples)

	## Turn the 2D array into a dataframe with multi indices
	df_density2D_IL = DataFrame(var_stacked.T,columns=m_index)

	## Save to csv
	df_density2D_IL.to_csv(os.path.join(outputdir,outputfile_density))










