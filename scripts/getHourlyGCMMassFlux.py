##############################################################################
#                                                                            #
#    Benjamin Fildier, Dec 2016.                                             #
#                                                                            #
#    This script extract hourly GCM-scale omega values along with the        #
#    pressure coordinates, compute the vertical integral of the mass flux    #
#    and stores it in a new variable MF.                                     #
#                                                                            #
##############################################################################


###--- Modules ---###
import glob
import sys,os
from netCDF4 import Dataset
from datetime import datetime
import numpy as np
import socket
from argparse import ArgumentParser

## Own functions
currentpath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,currentpath[:currentpath.rfind('/')+1]+'functions')

from importingData import *
from extremeScaling import *

###--- Functions ---###

def getInputfileNames(dates,compset,inputdir):

	dt_bnds = []
	for tstamp in dates:
		dt_info = [int(s) for s in tstamp.split('-') if s.isdigit()]
		dt_bnds.append(datetime(dt_info[0],dt_info[1],dt_info[2],dt_info[3]/3600))
	# Change date format to CMIP5 standard
	datetime1 = dt_bnds[0].isoformat().replace('-','').replace('T','').replace(':','')[:-2]
	datetime2 = dt_bnds[1].isoformat().replace('-','').replace('T','').replace(':','')[:-2]
	# Store correct files (between time boundaries, with correct compset name)
	ncfiles = []
	for file in glob.glob(os.path.join(inputdir,'*.nc')):
		filename = os.path.basename(file)
		dt_info = [int(s) for s in filename.split('.')[3].split('-') if s.isdigit()]
		dt = datetime(dt_info[0],dt_info[1],dt_info[2],dt_info[3]/3600)
		if compset in filename and ".cam.h0." in filename:
			if dt <= dt_bnds[1] and dt >= dt_bnds[0]:
				ncfiles.append(file)
	# Exit if no interesting files are found
	if len(ncfiles) == 0:
		print "No files matching request."
		sys.exit(0)
	ncfiles.sort()

	return ncfiles

def getInputDirectoryFx(inputdir):

	hostname = socket.gethostname()
	if hostname == "jollyjumper":
		dataroot = os.path.dirname(os.path.dirname(inputdir))
	elif "edison" in hostname:
		dataroot = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(inputdir))))
	inputdir_fx = os.path.join(dataroot,'preprocessed','allExperiments/fx')

	return inputdir_fx

def getOutputDirectory(inputdir,experiment,compset):

	case = "bf_"+compset+"_"+experiment
	currentpath = os.path.dirname(os.path.realpath(__file__))
	hostname = socket.gethostname()
	if hostname == "jollyjumper":
		outputdir = os.path.join(os.path.dirname(os.path.dirname(inputdir)),
			'preprocessed',case,'1hr')
	elif "edison" in hostname:
		outputdir = os.path.join(os.path.dirname(currentpath),'preprocessed',case,'1hr')

	return outputdir

def getDates(dates):
	# Get time boundaries
	dt_bnds = []
	for tstamp in dates:
		dt_info = [int(s) for s in tstamp.split('-') if s.isdigit()]
		dt_bnds.append(datetime(dt_info[0],dt_info[1],dt_info[2],dt_info[3]/3600))
	# Change date format to CMIP5 standard
	datetime1 = dt_bnds[0].isoformat().replace('-','').replace('T','').replace(':','')[:-2]
	datetime2 = dt_bnds[1].isoformat().replace('-','').replace('T','').replace(':','')[:-2]
	return datetime1, datetime2
	

def parseArguments():

	## Get arguments
	parser = ArgumentParser(description="Compute daily-mean integral of the vertical mass flux.")
	parser.add_argument('-e','--experiment',default="piControl",
		choices=['piControl','abrupt4xCO2'],help="Physical experiment or scenario.")
	parser.add_argument('-c','--compset',choices=['FSPCAMm_AMIP','FAMIPC5'],
		default="FSPCAMm_AMIP",help="CESM compset used.")
	parser.add_argument('-d','--dates',nargs=2,
		help="Time boundaries in format YYYY-MM-DD-sssss.")
	parser.add_argument('-i','--inputdir',help="Full path to input directory")
	args = parser.parse_args()

	## Get input directories
	inputdir_fx = getInputDirectoryFx(args.inputdir)
	print inputdir_fx

	## Format dates
	dates = getDates(args.dates)

	## Get outputdir and outputfile
	outputdir = getOutputDirectory(args.inputdir,args.experiment,args.compset)
	outputfile = "MF_1hr_CESM111-SPCAM20_"+args.experiment+"_r1i1p1_"+dates[0]+'-'+dates[1]+".nc"

	## Find input files
	ncfiles = getInputfileNames(args.dates,args.compset,args.inputdir)

	return args.experiment, args.compset, args.dates[0], args.dates[1], ncfiles, \
	args.inputdir, inputdir_fx, outputdir, outputfile


if __name__ == "__main__":


	##-- Get arguments --##

	experiment, case, datetime1, datetime2, ncfiles, inputdir, inputdir_fx, outputdir, \
	outputfile = parseArguments()

	varid = "OMEGA"

	lev_file = 'lev_fx_CESM111-SPCAM20_allExperiments_r0i0p0.nc'
	computeP = getPressureCoordinateFunction(os.path.join(inputdir_fx,lev_file))


	##-- Create output file --##	

	fh_source = Dataset(ncfiles[0],'r')
	lon_source = fh_source.variables['lon']
	lat_source = fh_source.variables['lat']
	vardims = fh_source.variables[varid].dimensions
	print vardims
	varshape = fh_source.variables[varid].shape
	rootgrp = Dataset(os.path.join(outputdir,outputfile),'w')

	#- Create dimensions and global attributes -#
	rootgrp.createDimension("lon",len(lon_source))
	nlat = len(lat_source)
	lat_slice = slice(nlat/3,2*nlat/3)
	rootgrp.createDimension("lat",nlat/3)
	rootgrp.createDimension("time",None)
	rootgrp.case = str(fh_source.case)
	rootgrp.description = "Hourly mass flux values from CESM111-SPCAM20 in the "+experiment+" scenario"
	
	#- Define variables -#
	lon = rootgrp.createVariable("lon","f8",("lon",))
	lon[:] = lon_source[:]
	lat = rootgrp.createVariable("lat","f8",("lat",))
	lat[:] = lat_source[lat_slice]
	time = rootgrp.createVariable("time","f8",("time",))
	date = rootgrp.createVariable("date","i4",("time",))
	datesec = rootgrp.createVariable("datesec","i4",("time",))
	var = rootgrp.createVariable('MF',"f4",('time','lat','lon'))

	#- Define variable attributes -#
	lon.long_name = str(lon_source.long_name)
	lat.long_name = str(lat_source.long_name)
	date.long_name = str(fh_source.variables['date'].long_name)
	datesec.long_name = str(fh_source.variables['datesec'].long_name)
	lon.units = str(lon_source.units)
	lat.units = str(lat_source.units)
	time.units = str(fh_source.variables['time'].units)
	time.calendar = str(fh_source.variables['time'].calendar)
	var.long_name = "Mean mass flux"
	var.units = str(fh_source.variables[varid].units)
	
	fh_source.close()

	##-- Read all files --##

	n = len(ncfiles)
	time_values = np.zeros((n,))
	date_values = np.zeros((n,))
	datesec_values = np.zeros((n,))
	var_values = np.zeros((n,nlat/3,varshape[-1]))
	
	# Start chronometer
	t0 = datetime.now()

	for i in range(n):

		ncfile = ncfiles[i]
		print "... Extracting "+ncfile.split('.')[3]+" data ..."
		# Open file
		fh = Dataset(ncfile)
		# Copy values
		time_values[i] = fh.variables['time'][:]
		date_values[i] = fh.variables['date'][:]
		datesec_values[i] = fh.variables['datesec'][:]
		# Get omega profile
		omega = fh.variables[varid][:,:,lat_slice,:]
		# get surface pressure
		pressurf = fh.variables["PS"][:,lat_slice,:]
		# Loop over all points
		massflux = np.zeros((1,nlat/3,varshape[-1]))
		for j in range(nlat/3):
			for k in range(varshape[-1]):
				pres = computeP(pressurf[0,j,k])
				massflux[0,j,k] = verticalPressureIntegral(pres,omega[0,:,j,k]) / verticalPressureIntegral(pres)
		var_values[i] = massflux
		# Close file
		fh.close()

	# Stop chronometer
	t1 = datetime.now()
	print
	print "------ It took", (t1-t0), "to process %d files."%n

	##-- Write to output --##

	time[:] = time_values[:]
	date[:] = date_values[:]
	datesec[:] = datesec_values[:]
	var[:] = var_values[:]

	rootgrp.close()

	sys.exit(0)
