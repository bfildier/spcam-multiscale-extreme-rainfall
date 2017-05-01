##############################################################################
#                                                                            #
#    Benjamin Fildier, Nov 2016.                                             #
#                                                                            #
#    Script extract the hourly values for a given varid, experiment name     #
#    from hourly files in a given directory, extract tropical values only    #
#    and write it to a single file.                                          #
#                                                                            #
#    Use         : varid defined on GCM grid (accept vertical coordinate)    #
#    Arguments   : varid experiment case datetime1 datetime2 inputdir        #
#    Output file : ${varid}_1hr_CESM111-SPCAM20_${experiment}_r1i1p1_...     #
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
## (varid, experiment, case, list of inputdir/*.nc files)
def parseArguments():
	args = sys.argv[1:]
	# Check number of arguments
	if len(args) != 6:
		print "Wrong number of arguments.\nRequired arguments: varid experiment case datetime1 datetime2 inputdir"
		sys.exit(1)
	# Get list of files in directory (3rd argument)
	varid, experiment, case, datetime1, datetime2 = args[:5]
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
	return varid, experiment, case, datetime1, datetime2, ncfiles


###--- Main program ---###

if __name__ == '__main__':

	##-- Get arguments --##

	varid, experiment, case, datetime1, datetime2, ncfiles = parseArguments()

	##-- Define input, output variables --##

	currentpath = os.path.dirname(os.path.realpath(__file__))
	inputdir = os.path.dirname(ncfiles[0])	
	hostname = socket.gethostname()
	if hostname == "jollyjumper":
		outputdir = os.path.join(os.path.dirname(os.path.dirname(inputdir)),
			'preprocessed',case,'1hr')
	elif "edison" in hostname:
		outputdir = os.path.join(os.path.dirname(currentpath),'preprocessed',case,'1hr')
	outputfile = varid.replace('_','-')+"_1hr_CESM111-SPCAM20_"+experiment+"_r1i1p1_"+datetime1+'-'+datetime2+".nc"

	##-- Create output file --##

	fh_source = Dataset(ncfiles[0],'r')
	lon_source = fh_source.variables['lon']
	lat_source = fh_source.variables['lat']
	vert = "lev" in fh_source.variables[varid].dimensions
	vardims = fh_source.variables[varid].dimensions
	varshape = fh_source.variables[varid].shape
	rootgrp = Dataset(os.path.join(outputdir,outputfile),'w')

	#- Create dimensions and global attributes -#
	rootgrp.createDimension("lon",len(lon_source))
	nlat = len(lat_source)
	lat_slice = slice(nlat/3,2*nlat/3)
	rootgrp.createDimension("lat",nlat/3)
	# If vertical coordinate exists
	if vert:
		rootgrp.createDimension("lev",len(fh_source.variables['lev']))
	rootgrp.createDimension("time",None)
	rootgrp.case = str(fh_source.case)
	rootgrp.description = "Hourly "+varid+" values from CESM111-SPCAM20 in the "+experiment+" scenario"
	
	#- Define variables -#
	lon = rootgrp.createVariable("lon","f8",("lon",))
	lon[:] = lon_source[:]
	lat = rootgrp.createVariable("lat","f8",("lat",))
	lat[:] = lat_source[lat_slice]
	# If vertical coordinate exists
	if vert:
		lev = rootgrp.createVariable("lev","f8",("lev",))
		lev[:] = fh_source.variables['lev'][:]
	time = rootgrp.createVariable("time","f8",("time",))
	date = rootgrp.createVariable("date","i4",("time",))
	datesec = rootgrp.createVariable("datesec","i4",("time",))
	var = rootgrp.createVariable(varid,"f4",vardims)

	#- Define variable attributes -#
	lon.long_name = str(lon_source.long_name)
	lat.long_name = str(lat_source.long_name)
	date.long_name = str(fh_source.variables['date'].long_name)
	datesec.long_name = str(fh_source.variables['datesec'].long_name)
	lon.units = str(lon_source.units)
	lat.units = str(lat_source.units)
	# If vertical coordinate exists
	if vert:
		lev.long_name = str(fh_source.variables['lev'].long_name)
		lev.units = str(fh_source.variables['lev'].units)
	time.units = str(fh_source.variables['time'].units)
	time.calendar = str(fh_source.variables['time'].calendar)
	var.long_name = str(fh_source.variables[varid].long_name)
	var.units = str(fh_source.variables[varid].units)
	
	fh_source.close()

	##-- Read all files --##

	n = len(ncfiles)
	time_values = np.zeros((n,))
	date_values = np.zeros((n,))
	datesec_values = np.zeros((n,))
	if vert:
		var_values = np.zeros((n,varshape[1],nlat/3,varshape[3]))
	else:
		var_values = np.zeros((n,nlat/3,varshape[2]))
	
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
		if vert:
			var_values[i] = fh.variables[varid][:,:,lat_slice,:]
		else:
			var_values[i] = fh.variables[varid][:,lat_slice,:]
		# Close file
		fh.close()

	# Stop chronometer
	t1 = datetime.now()
	# print "------------------------------------------------------------"
	print
	print "------ It took", (t1-t0), "to process %d files."%n
	# print "------------------------------------------------------------"

	##-- Write to output --##

	time[:] = time_values[:]
	date[:] = date_values[:]
	datesec[:] = datesec_values[:]
	var[:] = var_values[:]

	rootgrp.close()

	sys.exit(0)