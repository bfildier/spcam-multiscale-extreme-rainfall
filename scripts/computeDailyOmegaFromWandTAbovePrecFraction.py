##############################################################################
#                                                                            #
#    Benjamin Fildier, July 2017.                                            #
#                                                                            #
#    From (daily) convective-scale values of RHO-CRM-T-IXX and CRM_W_IXX,    #
#    compute the mean pressure velocity in the convective-scale events.      #
#    Call it OMEGA_CRM_WT_IXX.                                               #
#                                                                            #
##############################################################################


###--- Modules ---###
import sys,os,glob
from netCDF4 import Dataset
import numpy as np
import socket
from argparse import ArgumentParser
import datetime as dt
from datetime import datetime,timedelta

## Own functions
currentpath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,currentpath[:currentpath.rfind('/')+1]+'functions')

from importingData import *
from thermodynamics import airDensity,g


###--- Functions ---###

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

def parseArguments():

	parser = ArgumentParser(description="Extract GCM points corresponding to "+\
		"a given quantile of precipitation.")
	parser.add_argument('-f','--fraction',required=True,
		help="Threshold fraction alpha of rainfall amount to define rainfall event (%)."+\
			"If the integer entered is negative, all points are taken")
	parser.add_argument('-d','--dataroot',required=True,
		help="Input directory root.")
	parser.add_argument('-e','--experiment',default="piControl",
		choices=['piControl','abrupt4xCO2'],
		help="Physical experiment or scenario.")
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

	# Get case name
	case = "bf_%s_%s"%(args.compset,args.experiment)

	# Get input directories
	inputdir, inputdir_processed_day, inputdir_processed_1hr, inputdir_results, \
		inputdir_fx = getInputDirectories(args.dataroot,args.compset,args.experiment)

	# Store correct files (between time boundaries, with correct case name)
	ncfiles = []
	for file in glob.glob(inputdir+"/*.nc"):
		filename = os.path.basename(file)
		if case in filename and ".cam.h0." in filename:
			dt_info = [int(s) for s in filename.split('.')[3].split('-') if s.isdigit()]
			d = datetime(dt_info[0],dt_info[1],dt_info[2],dt_info[3]/3600)
			if d <= dt_bnds[1] and d >= dt_bnds[0]:
				ncfiles.append(file)
	# Exit if no interesting files are found
	if len(ncfiles) == 0:
		print "No files matching request."
		sys.exit(0)
	ncfiles.sort()

	return int(args.fraction), args.experiment, args.compset, args.dataroot, \
		datetime1[:8], datetime2[:8], ncfiles, args.dataroot


def defineOutputDirectory(dataroot,compset,experiment):

	hostname = socket.gethostname()
	case = "bf_%s_%s"%(compset,experiment)
	if hostname == "jollyjumper":
		outputdir = os.path.join(dataroot,'preprocessed',case,'day')
	elif "edison" in hostname or "cori" in hostname:
		outputdir = os.path.join(os.path.dirname(currentpath),'preprocessed',case,'day')

	return outputdir


###--- Main program ---###

if __name__ == "__main__":

	t0 = datetime.now()		# Start chronometer

	##-- Get arguments --##

	fraction, experiment, compset, case, datetime1, datetime2, ncfiles, \
		dataroot = parseArguments()

	##-- Get input directories --##

	inputdir, inputdir_processed_day, inputdir_processed_1hr, inputdir_results, \
		inputdir_fx = getInputDirectories(dataroot,compset,experiment)

	##-- Variables needed --##

	ref_id = "T"
	if fraction < 0:
		w_id = "W"
		rho_id = "RHO_MESO"
		omega_id = "OMEGA_MESO"
	else:
		indvarid = "I%d"%fraction
		w_id = 'CRM_W_%s'%indvarid
		rho_id = "RHO_CRM_T_%s"%indvarid
		omega_id = "OMEGA_CRM_WT_%s"%indvarid

	##-- Output metadata --##

	outputdir = defineOutputDirectory(dataroot,compset,experiment)
	outputfile = "%s_day_CESM111-SPCAM20_%s_r1i1p1_%s-%s.nc"%\
		(omega_id.replace('_','-'),experiment,datetime1,datetime2)

	##-- Create outputfile --##

	fh_source = Dataset(ncfiles[0],'r')
	lon_source = fh_source.variables['lon']
	lat_source = fh_source.variables['lat']
	lev_source = fh_source.variables['lev']
	omega_dims = fh_source.variables[ref_id].dimensions
	rootgrp = Dataset(os.path.join(outputdir,outputfile),'w')

		#- Create dimensions and global attributes -#

	rootgrp.createDimension("lon",len(lon_source))
	nlat = len(lat_source)
	lat_slice = slice(nlat/3,2*nlat/3)
	rootgrp.createDimension("lat",nlat/3)
	rootgrp.createDimension("lev",len(fh_source.variables['lev']))
	rootgrp.createDimension("time",None)
	rootgrp.case = str(fh_source.case)
	rootgrp.description = "Daily %s values from CESM111-SPCAM20 in the %s scenario"%\
		(omega_id,experiment)

		#- Define variables -#

	lon = rootgrp.createVariable("lon","f8",("lon",))
	lon[:] = lon_source[:]
	lat = rootgrp.createVariable("lat","f8",("lat",))
	lat[:] = lat_source[lat_slice]
	lev = rootgrp.createVariable("lev","f8",("lev",))
	lev[:] = fh_source.variables['lev'][:]
	time = rootgrp.createVariable("time","f8",("time",))
	date = rootgrp.createVariable("date","i4",("time",))
	datesec = rootgrp.createVariable("datesec","i4",("time",))
	omega = rootgrp.createVariable(omega_id,"f4",omega_dims)

	#- Define variable attributes -#
	lon.long_name = str(lon_source.long_name)
	lat.long_name = str(lat_source.long_name)
	lev.long_name = str(fh_source.variables['lev'].long_name)
	date.long_name = str(fh_source.variables['date'].long_name)
	datesec.long_name = str(fh_source.variables['datesec'].long_name)
	lon.units = str(lon_source.units)
	lat.units = str(lat_source.units)
	lev.units = str(fh_source.variables['lev'].units)
	time.units = str(fh_source.variables['time'].units)
	time.calendar = str(fh_source.variables['time'].calendar)
	omega.long_name = str("Pressure velocity")
	omega.units = str("Pa/s")
	
	fh_source.close()
	

	##-- Import variables --##

		#- Time variables -#

	processed_files = getInputfiles(ref_id,inputdir_processed_day,dates=(datetime1,datetime2))
	time_values = getVar("time",inputdir_processed_day,inputfiles=[processed_files[0],])
	date_values = getVar("date",inputdir_processed_day,inputfiles=[processed_files[0],])
	datesec_values = getVar("datesec",inputdir_processed_day,inputfiles=[processed_files[0],])

		#- Physical variables -#

	w_values = getVar(w_id,inputdir_processed_day,dates=(datetime1,datetime2))
	# Reverse w values along vertical dimension
	w_values = w_values[:,::-1,...]
	rho_values = getVar(rho_id,inputdir_processed_day,dates=(datetime1,datetime2))

	##-- Compute hourly rho profile at each GCM point --##

		#- Compute rho -#

	omega_values = -g*np.multiply(rho_values,w_values)

	##-- Write to output --##

	time[:] = time_values[:]
	date[:] = date_values[:]
	datesec[:] = datesec_values[:]
	omega[:] = omega_values[:]

	rootgrp.close()

	t1 = datetime.now()			# Stop chronometer
	print t1-t0

	sys.exit(0)



