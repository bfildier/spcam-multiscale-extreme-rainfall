


###--- Modules ---###

import glob
import sys,os
import string
import numpy as np
from datetime import date, time, datetime, timedelta
from netCDF4 import Dataset


###--- Functions ---###

## inputfiles is a list of inputfiles to open simultaneously.
## dates is a pair of date strings in the form 'YYYYMMDD'
## - if both inputfiles and dates are NA, use all files $inputdir/$varid_*
## - if inputfiles is given, use this list in priority
## - if inputfiles is not given but dates is given, use corresponding files
## between those dates
#
# 1. Get inputfiles in correct format
def getInputfiles(varid,inputdir,inputfiles=None,dates=None):
    ## Get the list of input files
    if inputfiles == None:
        if dates == None:    # Get all matching files
            inputfiles = glob.glob(os.path.join(inputdir,varid.replace('_','-')+"_*.nc"))
        else:    # Filter between dates
            dt_bnds = [datetime(int(d[:4]),int(d[4:6]),int(d[6:8])) for d in dates]
            inputfiles = []
            for file in glob.glob(os.path.join(inputdir,varid.replace('_','-')+"_*.nc")):
                filename = os.path.basename(file)
                dt_info = filename.split('.')[0].split('_')[-1].split('-')
                dt = [datetime(int(d[:4]),int(d[4:6]),int(d[6:8])) for d in dt_info]
                if dt[0] >= dt_bnds[0] and dt[1] <= dt_bnds[1]:
                    inputfiles.append(file)
    else:
        ## Append dirname to all files if necessary
        inputfiles = [os.path.join(inputdir,f) if inputdir not in f else f for f in inputfiles]
    inputfiles.sort()
    if len(inputfiles) == 0:
        print inputdir, varid, dates
    return inputfiles
#
# 2. Get values for processed data ($dataroot/preprocessed/$case/$freq/*)
def getVar(varid,inputdir,inputfiles=None,dates=None):
    if inputfiles is None:
        inputfiles = getInputfiles(varid,inputdir,inputfiles,dates)
    if len(inputfiles) == 0:
        print "Error %s: no matching input file."%varid
        return
    values_list = []
    if len(inputfiles) == 0:
        return np.array([])
    for file in inputfiles:
        fh = Dataset(file,'r')
        values_list.append(fh.variables[varid][:])
        fh.close()
    try:
        values_array = np.concatenate(values_list,axis=0)
        return values_array
    except ValueError:
        return values_list

## Get values from original hourly output
def getSimulationValues(varid,inputdir,dates):

    # Find valid inputfiles
    dt_bnds = [datetime(int(d[:4]),int(d[4:6]),int(d[6:8]),int(d[8:10]),int(d[10:12]))\
        for d in dates]
    inputfiles = []
    for file in glob.glob(os.path.join(inputdir,"*.nc")):
        filename = os.path.basename(file)
        dt_info = filename.split('.')[-2].replace('-','')
        dt = datetime(int(dt_info[:4]),int(dt_info[4:6]),int(dt_info[6:8]))+\
            timedelta(seconds=int(dt_info[8:13]))
        if dt >= dt_bnds[0] and dt <= dt_bnds[1]:
            inputfiles.append(file)
    inputfiles.sort()
    # Extract data
    if len(inputfiles) == 0:
        print "Error: no matching input file."
        return
    values_list = []
    for file in inputfiles:
        fh = Dataset(file,'r')
        if varid in fh.variables.keys():
            values_list.append(fh.variables[varid][:])
        fh.close()
    return np.concatenate(values_list,axis=0)

## Reads in file lev_fx_* with the required information
## Returns a function which takes a surface pressure value
## and returns a vector of pressure values
def getPressureCoordinateFunction(input_lev_file):
    fh = Dataset(input_lev_file,'r')
    hyam = fh.variables['hyam'][:]
    hybm = fh.variables['hybm'][:]
    P0 = fh.variables['P0'][:]
    fh.close()
    return lambda ps: (P0*hyam+ps*hybm)    # In hPa,mbar
