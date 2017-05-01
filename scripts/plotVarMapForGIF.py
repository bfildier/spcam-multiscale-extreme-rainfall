
###--- Modules ---###
from mpl_toolkits.basemap import Basemap
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import glob
import string
import datetime as dt
import sys,os
import socket

# Own modules
currentpath = os.path.dirname(os.path.realpath(__file__))


###--- Main program ---###
if __name__ == "__main__":

	plt.rcParams.update({'axes.labelsize': 'x-large',
						'axes.titlesize': 'x-large',
						'xtick.labelsize': 'large',
						'ytick.labelsize': 'large',
						'figure.subplot.top': '0.87',
						'figure.subplot.wspace': '0.5',
						'figure.subplot.hspace': '0.3'})

	##-- Get arguments --##

	if len(sys.argv) == 5:
		varid = sys.argv[1]
		startday = sys.argv[2]
		ndays = int(sys.argv[3])
		experiment = sys.argv[4]
	else:
		print "Error: wrong number of arguments."
		print "Required arguments: varid startday ndays experiment"
		sys.exit(-1)

	##-- Environment --##

	testing = False
	hostname = socket.gethostname()
	if hostname == "jollyjumper":
		if testing:
			inputdir_PI = '/Users/bfildier/Data/simulations/CESM111-SPCAM20/piControl/day/r1i1p1'
			inputdir_4xCO2 = '/Users/bfildier/Data/simulations/CESM111-SPCAM20/abrupt4xCO2/day/r1i1p1'	
		else:
			inputdir_PI = '/Volumes/LongTermStorage/Data/simulations/CESM111-SPCAM20/piControl/day/r1i1p1'
			inputdir_4xCO2 = '/Volumes/LongTermStorage/Data/simulations/CESM111-SPCAM20/abrupt4xCO2/day/r1i1p1'
	elif "edison" in hostname:	
		inputdir_PI = os.path.join(os.path.dirname(currentpath),'preprocessed/bf_FSPCAMm_AMIP_piControl/day')
		inputdir_4xCO2 = os.path.join(os.path.dirname(currentpath),'preprocessed/bf_FSPCAMm_AMIP_abrupt4xCO2/day')

	figdir = os.path.join(os.path.dirname(currentpath),'figures')

	if experiment == 'piControl':
		inputdir = inputdir_PI
	else:
		inputdir = inputdir_4xCO2

	##-- Draw maps --##

	for i in range(ndays):

		date_start = dt.datetime.strptime(startday,'%Y%m%d') + i*dt.timedelta(days=1)
		date_end = dt.datetime.strptime(startday,'%Y%m%d') + (i+1)*dt.timedelta(days=1)
		d1 = date_start.isoformat().split("T")[0].replace('-','')
		d2 = date_end.isoformat().split("T")[0].replace('-','')

		varfile = '%s_day_CESM111-SPCAM20_piControl_r1i1p1_%s-%s.nc'%(varid,d1,d2)
		print "plot", varfile

		fh = Dataset(os.path.join(inputdir,varfile),'r')
		temp = fh.variables[varid][0,:,:]

		lon1D = fh.variables['lon'][:]
		lat1D = fh.variables['lat'][:]
		lon2D, lat2D = np.meshgrid(lon1D,lat1D)
		fh.close()

		fig = plt.figure(figsize=(12,2))
		map = Basemap(projection='cyl',lat_0=0,lon_0=180,llcrnrlon=0,llcrnrlat=-30,urcrnrlon=357,urcrnrlat=30)
		map.contourf(lon2D,lat2D,temp,cmap=plt.cm.gist_ncar,levels=range(285,315))
		map.drawparallels(range(-90, 100, 30),labels=[1,0,0,1])
		map.drawmeridians(range(0,400,90),labels=[1,0,0,1])
		map.drawcoastlines()
		plt.title('Surface Temperature %s (K) - %s'%(varid,date_start.isoformat().split("T")[0][5:]))
		plt.colorbar(pad=0.02,fraction=0.085)

		plt.savefig(os.path.join(figdir,'maps','map_%s_%s.pdf'%(varid,date_start.isoformat().split("T")[0].replace('-',''))))
		plt.close()











