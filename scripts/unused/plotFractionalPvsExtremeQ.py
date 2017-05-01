



###--- Modules ---###
import os, sys
import matplotlib.pyplot as plt
import socket
import pickle


## Own functions
currentpath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,currentpath[:currentpath.rfind('/')+1]+'functions')

from importingData import *
from statistics import *
from extremeScaling import *
import argparse

def parseArguments():

	parser = argparse.ArgumentParser(description="Plot high quantiles vs. rank for a sequence of precipitation variables")
	parser.add_argument('pr_ids',nargs="+",
		help="Sequence of precipitation variables IDs to plot")
	parser.add_argument('-s','--subset',default=None,choices=['ocean','land'],
		help="Reduce the domain to a subset of points")
	parser.add_argument('-d','--directory',required=True,
		help="Input directory root")
	parser.add_argument('-c','--compset',choices=['FSPCAMm_AMIP','FAMIPC5'],
		default="FSPCAMm_AMIP",help="CESM compset used")
	args = parser.parse_args()

	return args.pr_ids, args.directory, args.compset, args.subset


def getInputDirectories(dataroot,compset):

	hostname = socket.gethostname()
	case_PI = "bf_%s_piControl"%compset
	case_4xCO2 = "bf_%s_abrupt4xCO2"%compset
	if hostname == "jollyjumper":
		inputdir_PI = os.path.join(dataroot,'preprocessed',case_PI,'day')
		inputdir_4xCO2 = os.path.join(dataroot,'preprocessed',case_4xCO2,'day')
	elif "edison" in hostname:	
		inputdir_PI = os.path.join(os.path.dirname(currentpath),
			'preprocessed',case_PI,'day')
		inputdir_4xCO2 = os.path.join(os.path.dirname(currentpath),
			'preprocessed',case_4xCO2,'day')
	inputdir_fx = os.path.join(os.path.dirname(os.path.dirname(inputdir_PI)),
		'allExperiments/fx')

	return inputdir_PI, inputdir_4xCO2, inputdir_fx

def getDeltaT(inputdir_PI,inputdir_4xCO2,subset):

	ts_PI = getVar('TS',inputdir_PI)
	ts_PI_m = computeMean(ts_PI,subset=subset)
	ts_4xCO2 = getVar('TS',inputdir_4xCO2)
	ts_4xCO2_m = computeMean(ts_4xCO2,subset=subset)
	dts = ts_4xCO2_m-ts_PI_m

	return dts

###--- Main program ---###
if __name__ == "__main__":

	##-- Get arguments --##

	pr_ids, dataroot, compset, subsetName = parseArguments()

	##-- Set environment --##

	inputdir_PI, inputdir_4xCO2, inputdir_fx = getInputDirectories(dataroot,compset)

	##-- Define subset if necessary --##

	if subsetName == "ocean":
		# Define subset of points over ocean
		landfile = 'landmask_fx_CESM111-SPCAM20_allExperiments_r0i0p0.nc'
		fh_landmask = Dataset(os.path.join(inputdir_fx,landfile))
		landmask = fh_landmask.variables['landmask']
		subset = np.logical_not(landmask)
		fh_landmask.close()
	elif subsetName == "land":
		# Define subset of points over land
		landfile = 'landmask_fx_CESM111-SPCAM20_allExperiments_r0i0p0.nc'
		fh_landmask = Dataset(os.path.join(inputdir_fx,landfile))
		landmask = fh_landmask.variables['landmask']
		subset = landmask
		fh_landmask.close()
	
	##-- Compute change in temperature --##
	
	dts = getDeltaT(inputdir_PI,inputdir_4xCO2,subsetName)
	print "Delta T =", dts

	##-- Import precipitation data --##

	# Initialize statistics
	pr_PI_c_IL_all = []
	pr_4xCO2_c_IL_all = []
	Nvars = len(pr_ids)

	for pr_id in pr_ids:

		##-- Import data --##

		pr_PI = getVar(pr_id,inputdir_PI)*86400*1000    # Convert to mm/day
		pr_4xCO2 = getVar(pr_id,inputdir_4xCO2)*86400*1000    # Convert to mm/day
		pr_units = "mm/day"

		##-- Compute statistics --##

		maxQ, Q_IL, pr_PI_c_IL, pr_PI_e_IL, pr_PI_c_L, pr_PI_e_L, pr_PI_Hd_L =\
			getStatistics(pr_PI,subset=subsetName)
		maxQ, Q_IL, pr_4xCO2_c_IL, pr_4xCO2_e_IL, pr_4xCO2_c_L, pr_4xCO2_e_L, pr_4xCO2_Hd_L =\
			getStatistics(pr_4xCO2,subset=subsetName)
		
		pr_PI_c_IL_all.append(pr_PI_c_IL)
		pr_4xCO2_c_IL_all.append(pr_4xCO2_c_IL)

	##-- Figure options --##

	plt.rcParams.update({'figure.subplot.left': '0.16',
						'figure.subplot.right': '0.94',
						'legend.fontsize':'small'})


	print "inputdir_PI    :", inputdir_PI
	print "inputdir_4xCO2 :", inputdir_4xCO2
	print "varids         :", pr_ids

	## Environment
	figdir = os.path.join(os.path.dirname(currentpath),'figures')

	## Colors
	colorfile = os.path.join(currentpath,'colorsAndTypes.pickle')
	with open(colorfile,'rb') as handle:
		col = pickle.load(handle)
		lt = pickle.load(handle)
		pal = pickle.load(handle)

	##-- Plot --#
	
	if subsetName is not None:
		figtitle = r"Fractional Change $\frac{\Delta P_q/P_q}{\Delta T_s} over %s$"%subsetName
		fignameroot = "fig_%s_fractionalPvsExtremeQ"%subsetName
	else:
		figtitle = r"Fractional Change $\frac{\Delta P_q/P_q}{\Delta T_s}$"
		fignameroot = "fig_fractionalPvsExtremeQ"


	fig, ((ax1)) = plt.subplots(ncols=1, nrows=1, figsize=(5.5,5))

	ax1.set_xscale('log')
	ax1.set_ylim((-5,20))
	# Trick to plot inverse logarithmic scale
	x = np.flipud(1./np.subtract(np.ones(Q_IL.size),Q_IL/100.))

	for (pr_PI_c_IL,pr_4xCO2_c_IL,pr_id) in zip(pr_PI_c_IL_all,pr_4xCO2_c_IL_all,pr_ids):
		print "plot %s"%pr_id
		# Plot curve
		frac_change_pr = 100.*(np.divide(pr_4xCO2_c_IL,pr_PI_c_IL)-np.ones(x.shape))/dts
		ax1.plot(x,frac_change_pr,label=pr_id,c=col[pr_id])
		# if Nvars > 1:
		# 	ax1.plot(x,frac_change_pr,label=pr_id)
		# else:
		# 	ax1.plot(x,frac_change_pr,c=col)

	ax1.invert_xaxis() # reverse x-axis
	labels = [item.get_text() for item in ax1.get_xticklabels()]
	n = ceil(log10(x.max()))
	N = len(labels)
	for i in range(1,N):
	    labels[i] = str(100*(1-10**(-n+i-1)))
	    if -n+i-1 == 0:
	        break
	ax1.set_xticklabels(labels)
	ax1.axhline(y=0,linestyle=':',c='k')
	ax1.set_title(figtitle)
	ax1.set_xlabel(r"Percentile $q$ (%)")
	ax1.set_ylabel(r"$\delta P /\Delta T_s$ (%/K)")
	ax1.legend(loc='upper left')

	if Nvars == 1:
		plt.savefig(os.path.join(figdir,'%s_%s.pdf'%(fignameroot,pr_id)))
		print '%s_%s.pdf'%(fignameroot,pr_id)
	else:
		plt.savefig(os.path.join(figdir,'%s_all.pdf'%fignameroot))
		print '%s_all.pdf'%fignameroot


