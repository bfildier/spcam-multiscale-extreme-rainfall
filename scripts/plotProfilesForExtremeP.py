###--- Modules ---###
import os, sys
import matplotlib.pyplot as plt
import socket

# Own modules
currentpath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,currentpath[:currentpath.rfind('/')+1]+'functions')

from importingData import *
from statistics import *
from thermodynamics import airDensity
from thermo_constants import *
g = gg


###--- Functions ---###

def computeMeanAndExtremeProfiles(pr_id,inputdir,quantile,computeP):
	omega_id = 'OMEGA'
	ts_id = 'TS'
	t_id = 'T'
	ps_id = 'PS'
	relhum_id = 'RELHUM'
	q_id = 'Q'
	# Extract all variables
	pr = getVar(pr_id,inputdir)
	omega = getVar(omega_id,inputdir)
	ts = getVar(ts_id,inputdir)
	t = getVar(t_id,inputdir)
	ps = getVar(ps_id,inputdir)
	relhum = getVar(relhum_id,inputdir)
	q = getVar(q_id,inputdir)
	# Compute means
	omega_m = computeMean(omega)
	ts_m = computeMean(ts)
	t_m = computeMean(t)
	ps_m = computeMean(ps)
	p_m = computeP(ps_m)
	relhum_m = computeMean(relhum)
	q_m = computeMean(q)
	rho_m = airDensity(t_m,p_m,q_m)
	w_m = np.divide(omega_m,-g*rho_m)
	profiles_m = (p_m,t_m,q_m,relhum_m,rho_m,omega_m,w_m)
	# Compute means for extreme quantile
	i_q = getIndicesOfExtremePercentile(pr,Q)
	omega_q = computeMeanAtTimeLatLonIndices(omega,i_q)
	ts_q = computeMeanAtTimeLatLonIndices(ts,i_q)
	t_q = computeMeanAtTimeLatLonIndices(t,i_q)
	ps_q = computeMeanAtTimeLatLonIndices(ps,i_q)
	p_q = computeP(ps_q)
	relhum_q = computeMeanAtTimeLatLonIndices(relhum,i_q)
	q_q = computeMeanAtTimeLatLonIndices(q,i_q)
	rho_q = airDensity(t_q,p_q,q_q)
	w_q = np.divide(omega_q,-g*rho_q)
	profiles_q = (p_q,t_q,q_q,relhum_q,rho_q,omega_q,w_q)
	# Return
	return profiles_m, profiles_q

def parseArguments():
	return sys.argv[1]

###--- Main program ---###
if __name__ == "__main__":

	plt.rcParams.update({'axes.labelsize': 'x-large',
                     'axes.titlesize': 'x-large',
                     'xtick.labelsize': 'large',
                     'ytick.labelsize': 'large',
                     'figure.subplot.top': '0.87',
                     'figure.subplot.wspace': '0.5',
                     'figure.subplot.hspace': '0.4',
                     'legend.fontsize':'medium'})

	##-- Arguments --##

	if len(sys.argv) == 4:
		pr_id = sys.argv[1]
		Q = float(sys.argv[2])
		dataroot = sys.argv[3]
	else:
		print "Error: wrong number of arguments."
		print "Required arguments: varid quantile"
		sys.exit(-1)

	##-- Import precipitation data --##

	hostname = socket.gethostname()
	case_PI = "bf_FSPCAMm_AMIP_piControl"
	case_4xCO2 = "bf_FSPCAMm_AMIP_abrupt4xCO2"
	if hostname == "jollyjumper":
		inputdir_PI = os.path.join(dataroot,'preprocessed',case_PI,'day')
		inputdir_4xCO2 = os.path.join(dataroot,'preprocessed',case_4xCO2,'day')
	elif "edison" in hostname:	
		inputdir_PI = os.path.join(os.path.dirname(currentpath),
			'preprocessed',case_PI,'day')
		inputdir_4xCO2 = os.path.join(os.path.dirname(currentpath),
			'preprocessed',case_4xCO2,'day')
	inputdir_fx = os.path.join(os.path.dirname(os.path.dirname(inputdir_PI)),'allExperiments/fx')

	print "inputdir_PI    :", inputdir_PI
	print "inputdir_4xCO2 :", inputdir_4xCO2
	print "varid          :", pr_id
	print "quantile       :", Q

	figdir = os.path.join(os.path.dirname(currentpath),'figures')

	##-- Compute profiles for PI control and 4xCO2
	
	# Vertical pressure coordinate
	lev_file = 'lev_fx_CESM111-SPCAM20_allExperiments_r0i0p0.nc'
	computeP = getPressureCoordinateFunction(os.path.join(inputdir_fx,lev_file))

	print "... Computing vertical profiles for piControl run ..."
	profiles_PI_m, profiles_PI_q = computeMeanAndExtremeProfiles(pr_id,inputdir_PI,Q,computeP)
	print "... Computing vertical profiles for abrupt4xCO2 run ..."
	profiles_4xCO2_m, profiles_4xCO2_q = computeMeanAndExtremeProfiles(pr_id,inputdir_4xCO2,Q,computeP)

	##-- Plotting --##

	print "... Plotting ..."
	xlabs = (r"Air temperature (K)",r"Specific humidity ()",r"Relative humidity (%)",
		r"Density (kg/m3)",r"Pressure velocity (Pa/s)",r"Vertical velocity (m/s)")
	tlabs = (r"$\overline{T}$ and $T\vert_{P=P_{%s}}$"%(str(Q)),
		r"$\overline{q}$ and $q\vert_{P=P_{%s}}$"%(str(Q)),
		r"$\overline{\mathcal{H}}$ and $\mathcal{H}\vert_{P=P_{%s}}$"%(str(Q)),
		r"$\overline{\rho}$ and $\rho\vert_{P=P_{%s}}$"%(str(Q)),
		r"$\overline{\omega}$ and $\omega\vert_{P=P_{%s}}$"%(str(Q)),
		r"$\overline{w}$ and $w\vert_{P=P_{%s}}$"%(str(Q)))

	fig, axs = plt.subplots(ncols=3, nrows=2, figsize=(15,10))
	fig.suptitle('Vertical profiles for the Tropics',fontsize='xx-large')

	for (ax,i) in zip(axs.ravel(),range(1,7)):
		ax.plot(profiles_PI_m[i],profiles_PI_m[0]/100.,'b--',label="Mean, PI")
		ax.plot(profiles_PI_q[i],profiles_PI_q[0]/100.,'b',label="For $P=P_{%s}$, PI"%(str(Q)))
		ax.plot(profiles_4xCO2_m[i],profiles_4xCO2_m[0]/100.,'r--',label="Mean, 4xCO2")
		ax.plot(profiles_4xCO2_q[i],profiles_4xCO2_q[0]/100.,'r',label="For $P=P_{%s}$, 4xCO2"%(str(Q)))
		ax.set_ylim(0,1000)
		ax.invert_yaxis()
		ax.set_title(tlabs[i-1])
		ax.set_xlabel(xlabs[i-1])
		ax.set_ylabel("Pressure (hPa)")
		if i == 1:
			ax.legend(loc='lower left')

	plt.savefig(os.path.join(figdir,'fig_profiles_PI_4xCO2_mean_extreme_%s_%s.pdf'%(str(Q).replace('.','-'),pr_id)))

	


