#################################
# Initialize physical constants #
#################################

from math import pi

# Earth
g = 9.81                     # Gravity; m2/s
r_earth = 6.37e6             # Earth's radius; m
Omega_earth = 2*pi/86400     # Earth's angular rotation; rad/sec

# Radiation
sigma_stefan = 5.67e-8       # Stefan-Boltzmann constant; W/m2/K4

# Moist thermodynamics
R_d = 287.04                 # Gas constant for dry air; J/kg/K
R_v = 461.50                 # Gas constant for water vapor; J/kg/K
c_p = 1005.7                 # Heat capacity at constant pressure of dry air; J/kg/K
c_pv = 1870                  # Heat capacity at constant pressure of water vapor; J/kg/K
c_pw = 4190                  # Heat capacity at constant pressure of liquis water; J/kg/K
L = 2.5e6                    # Latent heat of vaporisation at 0C; J/kg
# eps = 0.622                  # molwt_water/"molwt"_dryair=18.0/28.9
# delta = 0.608                # 1/eps-1, coefficient you get after linearizing (1+r/eps)/(1+r)