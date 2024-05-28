#%% IMPORTS
from matplotlib import pyplot as plt
import numpy as np
import json
#%% FILE HANDLING
# Opening JSON files
with open('data.json') as f:
    data = json.load(f)
    
# Organizing sets of data
transistor_details = {key: val for detail in data['transistor_details'] for key, val in detail.items()}
fin_constants = {key: val for detail in data['fin_constants'] for key, val in detail.items()}
system_details = {key: val for detail in data['system_details'] for key, val in detail.items()}
#%% DEFINING VARIABLES
N = np.linspace(1, 100, 1E-2)
D = np.linsapce(1E-5, 100, 1E-3)
L = np.linsapce(1E-5, 1E-3, 1E-7)

# DIMENSIONS OF THE TRANSISTOR SURFACE
pi = np.pi # defined it this way so it's easier to type for equations
conversion = 1E-3 # conversion from millimeters to meters
H = transistor_details['height'] * conversion # (m) height of the transistor
W = transistor_details['width'] * conversion # (m) width of the transistor

# SYSTEM PARAMETERS
density = fin_constants['density'] # density of the fin in kg/m^3
overdesign_constraint = 0.2 # percent overdesign taken into account
q = system_details['total_heat_loss'] * overdesign_constraint # (W) heat loss of the system 
Tw = system_details['surface_temp'] # (Celsius) top surface temperature that must be maintained
T_inf = system_details['environment_temp'] # (Celsius) surrounding environment temperature
dt = Tw - T_inf # (Celsius) temperature difference between surface and environment
h = system_details['film_heat_transfer_coefficient'] # (W/m^2*C) film heat transfer coefficient

# THERMAL CONDUCTIVITY - interpolating the values
slope = (fin_constants['thermal_conductivity@100'] - fin_constants['thermal_conductivity@0']) / (fin_constants['thermal_temp100'] - fin_constants['thermal_temp0'])
k = slope * (Tw - fin_constants['thermal_temp0']) + fin_constants['thermal_conductivity@0'] # (W/m*C) FINAL THERMAL CONDUCTIVITY AT 85 CELSIUS
#%% FUNCTIONS
def fin_effectiveness(N, D, L): # how effective the fin actually is
    P = D * np.pi # perimeter of the fin
    A_f = P * L # surface area of the fin
    
    m = np.sqrt((h * P) / (k * A_f)) # m-value
    mL = m * L
    if mL > 700:
        tanh_mL = 1
    elif mL < -700:
        tanh_mL = -1
    else:
        tanh_mL = np.tanh(mL)
        
    return tanh_mL / m
#%% PLOTTING
effectiveness = fin_effectiveness(N, D, L)
plt.plot(N, effectiveness)
plt.plot(D, effectiveness)
plt.plot(L, effectiveness)
plt.show() # shows the figure