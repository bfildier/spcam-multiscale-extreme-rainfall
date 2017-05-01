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
	parser.add_argument('-q','--qRanks',nargs=2,required=True,help='quantile')
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
		help="Boundary datetime values in YYYY-MM-DD-SSSSS format.")
	args = parser.parse_args()

	# Get time boundaries
	dt_bnds = []
	for tstamp in datetimes:
		dt_info = [int(s) for s in tstamp.split('-') if s.isdigit()]
		dt_bnds.append(datetime(dt_info[0],dt_info[1],dt_info[2],dt_info[3]/3600))
	# Change date format to CMIP5 standard
	datetime1 = dt_bnds[0].isoformat().replace('-','').replace('T','').replace(':','')[:-2]
	datetime2 = dt_bnds[1].isoformat().replace('-','').replace('T','').replace(':','')[:-2]

	return args.pr_id, args.qRanks, args.directory, args.compset, args.experiment,
		args.subset, (datetime1,datetime2)

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

def getPercentileFile(referencepath,compset,experiment,subset):

	precdir = os.path.join(os.path.dirname(referencepath),'results')
	precfile = "pr_density_IL_%s_%s_%s.csv"%(compset,experiment,subset)

	return precdir,precfile

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


###---- Main program ---###

if __name__ == "__main__":

	##-- Get arguments --##

	pr_id, qRanks, dataroot, compset, experiment, subset, datetimes =  parseArguments()
	newvarid = "J%s-%s"%(string.join(qRanks[0].split('.')),pr_id)

	##-- Set environment --##

	inputdir, inputdir_fx = getInputDirectories(dataroot,compset,experiment)
	perdir, percfile = getPercentileFile(currentpath,compset,experiment,subset)
	outputdir = inputdir
	outputfile = "%s_day_CESM111-SPCAM20_%s_r1i1p1_%s-%s.nc"

	print pr_id, qRanks, dataroot, compset, experiment, subset, datetimes
	print outputfile

	##-- Get precipitation range values


	##-- Get daily precipitation values

	##-- Define subset if necessary --##

	##-- Create outputfile --##

	# binary variable with same dimensions

	##-- Compute mask of correct indices, intersect with subset and store to output --##

















