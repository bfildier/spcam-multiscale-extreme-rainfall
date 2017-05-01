#########################################################################
#                                                                       #
#    Courtesy of Jake Seeley, Nov 2016.                                 #
#                                                                       #
#########################################################################

##  Dry air                                                            ##
c_vd             = 719.             # J/kg/K                           ##
R_d              = 287.04           # J/kg/K                           ##
c_pd             = c_vd + R_d       # J/kg/K                           ##
                                                                       ##
##  Water vapor                                                        ##
c_vv             = 1418.            # J/kg/K                           ##
R_v              = 461.4            # J/kg/K                           ##
c_pv             = c_vv + R_v       # J/kg/K                           ##
                                                                       ##
##  Liquid water                                                       ##
c_vl             = 4216.            # J/kg/K                           ##
                                                                       ##
##  Solid water                                                        ##
c_vs             = 2106.            # J/kg/K                           ##    
                                                                       ##
##  Reference temperatures and pressures                               ##
T_0              = 273.16           # K                                ## 
p_0              = 1.e5             # Pa                               ##
e_0              = 611.65           # Pa                               ##
                                                                       ##
##  Energies, enthalpies, entropies                                    ##
L_v              = 2500000.         # J/kg                             ##
E_0v             = 2374000.         # J/kg                             ##
E_0s             = 333700.          # J/kg                             ##
s_0v             = E_0v/T_0 + R_v   # J/kg/K                           ##
s_0s             = E_0s/T_0         # J/kg/K                           ##
                                                                       ##
##  other                                                              ##
gg               = 9.81             # m/s^2                            ##
eps              = R_v/R_d - 1.     # Unitless                         ##
                                                                       ##
                                                                       ##
#########################################################################