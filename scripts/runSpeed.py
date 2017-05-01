###########################################################################
#                                                                         #
#    Benjamin Fildier, Nov 2016.                                          #
#    Plots the run cost as a function of number of nodes on lawrencium    #
#                                                                         #
###########################################################################

import os, sys
import matplotlib.pyplot as plt
currentpath = os.path.dirname(os.path.realpath(__file__))

fig = plt.figure(figsize=(5,5))
PEs = (5,10,15,20,25)
cost = (24521, 25283.45, 26053.05, 27740.17, 29431.69)
plt.plot(PEs,cost,'p')
plt.xlim(0,30)
plt.xlabel('# nodes')
plt.ylabel('Cost (PE.hrs/m-year)')
plt.title('on Lawrencium')
fig.savefig(os.path.join(os.path.dirname(currentpath),
	'figures/cost_SPCAM_on_lawrencium.pdf'),bbox_inches='tight')