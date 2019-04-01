##############################################################################
#                                                                            #
#    Benjamin Fildier, Nov 2016.                                             #
#                                                                            #
#    This script extracts, from SPCAM hourly outputs, values for var         #
#    CRM_PREC, for a given experiment name and input directory and time      #
#    boundaries. It computes the CRM locations where rainfall contributes    #
#    a fraction XX% of its largest intensities.                              #
#    Output variable: IXX (e.g. I90 for the most intense 90% rain)           #
#                                                                            #
#    Use         : CRM_PREC (time,crm_ny,crm_nz,lat,lon)                     #
#    Arguments   : fraction experiment case datetime1 datetime2 inputdir     #
#    Output file : I${fraction}_1hr_CESM111-SPCAM20_${experiment}_r1i1p1_... #
#                                          ...${datetime1}-${datetime2}.nc   #
#                                                                            #
##############################################################################


###--- Modules ---###
import glob
import sys,os
from netCDF4 import Dataset
from datetime import datetime
import numpy as np
import socket


###--- Functions ---###

## Check we get correct arguments and return the tuple:
## (fraction, experiment, case, datetime1, datetime2, list of inputdir/*.nc files)
def parseArguments():
	args = sys.argv[1:]
	# Check number of arguments
	if len(args) != 6:
		print "Wrong number of arguments.\nRequired arguments: fraction experiment case datetime1 datetime2 inputdir"
		sys.exit(1)
	# Get list of files in directory (3rd argument)
	fraction, experiment, case, datetime1, datetime2 = args[:5]
	# Get time boundaries
	dt_bnds = []
	for tstamp in [datetime1,datetime2]:
		dt_info = [int(s) for s in tstamp.split('-') if s.isdigit()]
		dt_bnds.append(datetime(dt_info[0],dt_info[1],dt_info[2],dt_info[3]/3600))
	# Change date format to CMIP5 standard
	datetime1 = dt_bnds[0].isoformat().replace('-','').replace('T','').replace(':','')[:-2]
	datetime2 = dt_bnds[1].isoformat().replace('-','').replace('T','').replace(':','')[:-2]
	# Store correct files (between time boundaries, with correct case name)
	ncfiles = []
	for file in glob.glob(args[-1]+'/*.nc'):
		filename = os.path.basename(file)
		dt_info = [int(s) for s in filename.split('.')[3].split('-') if s.isdigit()]
		dt = datetime(dt_info[0],dt_info[1],dt_info[2],dt_info[3]/3600)
		if case in filename and ".cam.h0." in filename:
			if dt <= dt_bnds[1] and dt >= dt_bnds[0]:
				ncfiles.append(file)
	# Exit if no interesting files are found
	if len(ncfiles) == 0:
		print "No files matching request."
		sys.exit(0)
	ncfiles.sort()
	return int(fraction), experiment, case, datetime1, datetime2, ncfiles

## Get boolean ndarray of locations contributing to the largest $fraction
## fraction of values
def indicesOfValuesAboveFraction(values,fraction=0.9):
	if fraction == 1:
		values_min = 0.
	else:
		values_sorted = np.flipud(np.sort(values.flatten()))
		values_cum = np.cumsum(values_sorted)
		i_max = np.argmax(values_cum > fraction*values_sorted.sum())
		values_min = values_sorted[i_max]
	return np.reshape(values > values_min, values.shape)


###--- Main program ---###

if __name__ == '__main__':

	##-- Get arguments --##

	fraction, experiment, case, datetime1, datetime2, ncfiles = parseArguments()
	varid = 'CRM_PREC'
	indvarid = "I%d"%fraction

	##-- Define input, output variables --##

	currentpath = os.path.dirname(os.path.realpath(__file__))
	inputdir = os.path.dirname(ncfiles[0])
	hostname = socket.gethostname()
	if hostname == "jollyjumper":
		outputdir = os.path.join(os.path.dirname(os.path.dirname(inputdir)),
			'preprocessed',case,'1hr')
	elif "edison" in hostname or "cori" in hostname:
		outputdir = os.path.join(os.path.dirname(currentpath),'preprocessed',case,'1hr')
	outputfile = indvarid+"_1hr_CESM111-SPCAM20_"+experiment+"_r1i1p1_"+datetime1+'-'+datetime2+".nc"

	##-- Create output file --##

	fh_source = Dataset(ncfiles[0],'r')
	lon_source = fh_source.variables['lon']
	lat_source = fh_source.variables['lat']
	vardims = fh_source.variables[varid].dimensions
	varshape = fh_source.variables[varid].shape
	rootgrp = Dataset(os.path.join(outputdir,outputfile),'w')

	#- Create dimensions and global attributes -#
	nlon = len(lon_source)
	rootgrp.createDimension("lon",nlon)
	nlat = len(lat_source)
	lat_slice = slice(nlat/3,2*nlat/3)
	rootgrp.createDimension("lat",nlat/3)
	rootgrp.createDimension("crm_nx",varshape[2])
	rootgrp.createDimension("crm_ny",varshape[1])
	rootgrp.createDimension("time",None)
	rootgrp.case = str(fh_source.case)
	rootgrp.description = "CRM Points contributing to the most intense "+\
	str(fraction)+"%% CRM_PREC values in the CESM111-SPCAM20 "+experiment+\
	" AMIP experiment."

	#- Define variables -#
	lon = rootgrp.createVariable("lon","f8",("lon",))
	lon[:] = lon_source[:]
	lat = rootgrp.createVariable("lat","f8",("lat",))
	lat[:] = lat_source[lat_slice]
	time = rootgrp.createVariable("time","f8",("time",))
	date = rootgrp.createVariable("date","i4",("time",))
	datesec = rootgrp.createVariable("datesec","i4",("time",))
	indvar = rootgrp.createVariable(indvarid,"u1",vardims)

	#- Define variable attributes -#
	lon.long_name = str(lon_source.long_name)
	lat.long_name = str(lat_source.long_name)
	date.long_name = str(fh_source.variables['date'].long_name)
	datesec.long_name = str(fh_source.variables['datesec'].long_name)
	lon.units = str(lon_source.units)
	lat.units = str(lat_source.units)
	time.units = str(fh_source.variables['time'].units)
	time.calendar = str(fh_source.variables['time'].calendar)
	indvar.long_name = "Location of the largest %d%% CRM_PREC values"%fraction
	indvar.units = ""

	fh_source.close()

	##-- Read all files --##

	n = len(ncfiles)
	time_values = np.zeros((n,))
	date_values = np.zeros((n,))
	datesec_values = np.zeros((n,))
	var_values = np.zeros((n,varshape[1],varshape[2],nlat/3,nlon))

	# Start chronometer
	t0 = datetime.now()

	# Read in all CRM_PREC values
	for i in range(n):

		ncfile = ncfiles[i]
		print "... Extracting "+ncfile.split('.')[3]+" data ..."
		# Open file
		fh = Dataset(ncfile)
		# Copy values
		time_values[i] = fh.variables['time'][:]
		date_values[i] = fh.variables['date'][:]
		datesec_values[i] = fh.variables['datesec'][:]
		var_values[i] = fh.variables[varid][:,:,:,lat_slice,:]
		# Close file
		fh.close()

	# Get indices of most intense XX% CRM_PREC
	indvar_values = np.zeros((n,varshape[1],varshape[2],nlat/3,nlon))
	for ilat in range(nlat/3):
		for ilon in range(nlon):
			values = var_values[:,:,:,ilat,ilon]
			f = float(fraction)/100.
			indvar_values[:,:,:,ilat,ilon] = indicesOfValuesAboveFraction(values,f)

	# Stop chronometer
	t1 = datetime.now()

	print
	print "------ It took", (t1-t0), "to process %d files."%n

	##-- Write to output --##

	time[:] = time_values[:].mean()
	date[:] = date_values[0]
	datesec[:] = datesec_values[:].mean()
	indvar[:] = indvar_values[:]

	rootgrp.close()

	sys.exit(0)
