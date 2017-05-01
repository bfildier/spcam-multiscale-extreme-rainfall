


## Modules
from matplotlib import colors
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
import pickle
import numpy as np

## Save file
savefile = 'colorsAndTypes.pickle'


##--  Define single colors  --##

col = {}

## For pr type PRECT
col['PRECT'] = 'k'

## For compsets
col['FSPCAMm_AMIP'] = 'k'
col['FAMIPC5'] = 'b'

## For subsets
col['tropics'] = 'k'
col['ocean'] = 'darkcyan'
col['land'] = 'saddlebrown'

##--  Define line types  --##

lt = {}

# For experiments
lt['piControl'] = '--'
lt['abrupt4xCO2'] = '-'

# For fracChange vs. scaling
lt['simulated'] = '-'
lt['scaling'] = ':'


##-- Define line widths --##
lw = {}

# For fracChange vs. scaling
lw['simulated'] = 1
lw['scaling'] = 2


##-- Define dashes --##
dashes = {}

# For contributions to fractional changes in pr
dashes['mass'] = [2, 8]
dashes['thermo'] = [8, 4]
dashes['dyn_mean'] = [8, 4, 8, 4, 2, 4, 2, 4]
dashes['dyn_mean_th'] = [8, 4, 8, 4, 4, 4, 4, 4]
dashes['dyn_mean_dyn'] = [8, 4, 6, 4, 4, 4, 2, 4]
dashes['dyn_profile'] = [8, 4, 2, 4, 2, 4, 2, 4, 2, 4]
dashes['pressure'] = [8, 8, 2, 8]
dashes['non_linear'] = [2, 4]
dashes['all'] = []

##-- Define markers --##
mark = {}

# For experiments
mark['piControl'] = 'o'
mark['abrupt4xCO2'] = "s"


##--  Define color palettes  --##
cmaps = {}
pal = {}

## CRM_PREC_IXX for 5 values (90,75,50,25,10)
pal_crm_prec = sns.cubehelix_palette(5, start=2, rot=-1.5, dark=0.6, light=0.3,
	gamma=1, hue=3, reverse=False)
pal['CRM_PREC'] = {}
crm_prec_names = ['CRM_PREC_I90','CRM_PREC_I75','CRM_PREC_I50','CRM_PREC_I25','CRM_PREC_I10']
for i in range(len(pal_crm_prec)):
	pal['CRM_PREC'][crm_prec_names[i]] = pal_crm_prec[i]
	col[crm_prec_names[i]] = pal_crm_prec[i]

## Percentile ranks for 7 values (90,96.84,99,99.68,99.9,99.97,99.99)
cm = plt.get_cmap('gnuplot_r')
cmaps["Q_IL"] = 'gnuplot_r'
# Q_names = ["90","96.84","99","99.684","99.9","99.968","99.99","99.997"]
Q_names = ['90.0', '92.0567', '93.6904', '94.9881', '96.0189', '96.8377',
       '97.4881', '98.0047', '98.4151', '98.7411', '99.0', '99.2057',
       '99.369', '99.4988', '99.6019', '99.6838', '99.7488', '99.8005',
       '99.8415', '99.8741', '99.9', '99.9206', '99.9369', '99.9499',
       '99.9602', '99.9684', '99.9749', '99.98', '99.9842', '99.9874',
       '99.99', '99.9921', '99.9937', '99.995', '99.996', '99.9968',
       '99.9975', '99.998', '99.9984', '99.9987','99.999']
pal_Qs = [cm(i) for i in np.linspace(0,1,len(Q_names))]
pal['Q_IL'] = {}
for i in range(len(pal_Qs)):
	pal['Q_IL'][Q_names[i]] = pal_Qs[i]
	col[Q_names[i]] = pal_Qs[i]

##--  Save all  --##
print "Saving all variables to %s"%savefile
vars = [col,lt,pal,lw,dashes,mark,cmaps]
with open(savefile,'wb') as handle:
	for var in vars:
		pickle.dump(var,handle)

