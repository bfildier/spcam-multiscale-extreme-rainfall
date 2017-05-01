################################################################################
#                                                                              #
#    Benjamin Fildier, Nov 2016.                                               #
#                                                                              #
#    Compute the mean profile for 3D variables on the GCM grid, for quantiles  #
#    of given precipitation variables: for all varids, pr_ids and Qs given in  #
#    input, and save them in a data frame.                                     #
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

	parser = ArgumentParser(description="Compute mean vertical profiles for a sequence of 3D variables\n"+
		"at locations of given percentiles of precipitation variables.")
	parser.add_argument('varids',nargs="+",
		help="Sequence of 3D variables defined on GCM grid.")
	parser.add_argument('-p','--pr_ids',nargs="+",required=True,
		help="Sequence of precipitation variables used to compute location of quantiles.")
	parser.add_argument('-q','--qRanks',nargs="+",required=True,
		help="Sequence of quantile ranks.")
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
	args = parser.parse_args()

	return args.varids, args.pr_ids, args.qRanks, args.directory, args.compset, args.experiment, args.subset

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

##-- To turn into dataframe --##

def stackInVaridAndPridDimensions(var):

	if np.array(var).shape[0] == 1:
		var_stacked_once = np.array(var[0])
	else:
		var_stacked_once = np.vstack(var)
	if np.array(var_stacked_once).shape[0] == 1:
		var_stacked = var_stacked_once[0]
	else:
		var_stacked = np.hstack(var_stacked_once)

	return var_stacked

def processNones(var_stacked,removeEmptyProfiles=True):

	nlev = len(next(v for v in var_stacked if v is not None))
	ind_to_remove = []
	for i in range(len(var_stacked)):
		if var_stacked[i] is None:
			if removeEmptyProfiles:    # then store relevant indices
				ind_to_remove.append(i)
			else:                      # then turn into a list of Nones
				var_stacked[i] = [None]*nlev

	return var_stacked, ind_to_remove


###--- Main program ---###

if __name__ == "__main__":

	##-- Get arguments --##

	varids, pr_ids, qRanks, dataroot, compset, experiment, subset = parseArguments()

	##-- Set environment --##

	inputdir, inputdir_fx = getInputDirectories(dataroot,compset,experiment)
	outputdir = os.path.join(os.path.dirname(currentpath),'results')
	outputfile = "var3D_meanProfileAtQ_%s_%s_%s.csv"%(compset,experiment,subset)

	##-- Define subset if necessary --##

	subset_pts = getSubset(subset,inputdir_fx,experiment=experiment)

	##-- Compute means --##

	var_meanProfileAtQ = []
	varids_to_remove = []
	pr_ids_to_remove = []

	for varid in varids:

		print "  Target variable:", varid

		##-- Import target variable data --##

		var = getVar(varid,inputdir)
		if var is None:
			print varid, "ignored."
			varids_to_remove.append(varid)
			continue
		elif len(var.shape) == 3:
			print varid, "is 2D, ignored."
			varids_to_remove.append(varid)
			continue
		else:
			i = np.argmax(np.array(varids) == varid) - len(varids_to_remove)

		var_meanProfileAtQ.append([])
		pr_ids_to_remove = []

		for pr_id in pr_ids:

			print "    Reference precipitation variable:", pr_id, ";",

			##-- Import reference precipitation data --##

			pr = getVar(pr_id,inputdir)
			if pr is None:
				print pr_id, "ignored."
				pr_ids_to_remove.append(pr_id)
				continue
			else:
				j = np.argmax(np.array(pr_ids) == pr_id) - len(pr_ids_to_remove)

			var_meanProfileAtQ[i].append([])

			for q in qRanks:

				Q = float(q)

				print q,

				##-- Find quantile locations --##

				i_q = getIndicesOfExtremePercentile(pr,Q,subset=subset_pts,n_unit_min=3)

				##-- Compute means at quantiles --##

				profileAtQ = computeMeanAtTimeLatLonIndices(var,i_q)
				if profileAtQ is not None:
					profileAtQ = profileAtQ.tolist()
				var_meanProfileAtQ[i][j].append(profileAtQ)

			print
		
		print

	## Remove incorrect variables
	for varid in varids_to_remove:
		varids.remove(varid)
	for pr_id in pr_ids_to_remove:
		pr_ids.remove(pr_id)

	print "Remaining target variables:", varids
	print "Remaining precipitation variables:", pr_ids

	##-- Turn results into dataframe --##

	## Stack twice to get a 1D array of lists
	var_stacked = stackInVaridAndPridDimensions(var_meanProfileAtQ)

	## Look for None values in the sequence, 
	## then remove them or turn then into a list of Nones with same length
	var_stacked, ind_to_remove = processNones(var_stacked,removeEmptyProfiles=True)

	## Remove empty profiles (sequences of Nones)
	var_stacked = var_stacked.tolist()
	for i in ind_to_remove[::-1]:
		del var_stacked[i]

	## Turn it into a 2D array
	var_stacked_2D = np.vstack(var_stacked)

	## Define multidimensional indices for dataframe
	vars = [[varid]*len(pr_ids)*len(qRanks) for varid in varids]
	vars = np.hstack(vars)
	prs = [[pr_id]*len(qRanks) for pr_id in pr_ids]
	prs = np.hstack(prs)
	prs = np.hstack([prs]*len(varids))
	Qs = qRanks*len(varids)*len(pr_ids)

	## Remove multi-indices of empty profiles 
	ind_tuples = zip(vars,prs,Qs)
	for i in ind_to_remove[::-1]:
		del ind_tuples[i]
	m_index = MultiIndex.from_tuples(ind_tuples)

	## Turn the 2D array into a dataframe with multi indices
	df_meanAtQ = DataFrame(var_stacked_2D.T,columns=m_index)

	##-- Save to CSV file --##

	df_meanAtQ.to_csv(os.path.join(outputdir,outputfile))





