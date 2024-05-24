#%% IMPORTS
import numpy as np
from scipy.optimize import minimize
import math
import json

#%% FILE HANDLING
# Opening JSON files
with open('data.json') as f:
    data = json.load(f)
    
# Organizing sets of data
transistor_details_RAW = data['transistor_details']
fin_constants_RAW = data['fin_constants']
system_details_RAW = data['system_details']

# Initializing the dictionaries
transistor_details = {}
fin_constants = {}
system_details = {}

# Details of the data sets for use
for detail in transistor_details_RAW:
    transistor_details.update(detail)
for detail in fin_constants_RAW:
    fin_constants.update(detail)
for detail in system_details_RAW:
    system_details.update(detail)

#%% DEFINING KNOWN SYSTEM VARIABLES
# DIMENSIONS OF THE TRANSISTOR SURFACE
conversion = 1E-3 # conversion from milimeters to meters
H = transistor_details['height']*conversion # (m) height of the transistor
W = transistor_details['width']*conversion # (m) width of the transistor

# SYSTEM PARAMETERS
density = fin_constants['density'] # density of the fin in kg/m3
overdesign_constraint = 0.2 # percent overdesign taken into account
q = system_details['total_heat_loss']*overdesign_constraint # (W) heat loss of the system 
Tw = system_details['surface_temp'] # (Celsius) top surface temperature that must be maintained
T_inf = system_details['environment_temp'] # (Celsius) surroudning environment temperature

# HEAT TRANSFER COEFFICIENT
h = system_details['film_heat_transfer_coefficient'] # (W/m^2*C) film heat transfer coefficient

# THERMCAL CONDUCTIVITY
# Interpolating for the proepr k-value of aluminum at 85 celsius
slope = (fin_constants['thermal_conductivity@100']-fin_constants['thermal_conductivity@0']) / (fin_constants['thermal_temp100'] - fin_constants['thermal_temp0'])
k = slope*(Tw-fin_constants['thermal_temp0']) + fin_constants['thermal_conductivity@0'] # (W/m*C) FINAL THERMAL CONDUCTIVY AT 85 CELSIUS

#%% INITIALING UNKNOWN PARAMETERS
diameter = 1E-10 # (m) diameter of the fin
length = 1E-10 # (m) length of the fin
num_of_fins = 1 # number of fins total

# To make things easier, grouping variables together to make equation less long
pi = math.pi # just to shorten it so we don't have to type out the 'math' part
m = math.sqrt((4*h)/(diameter*k)) # m-parameter
dt = Tw - T_inf # (Celsius) temperature difference between surface and environment

#%% MIN-MAX CALCULATION
# The goal is to minimize the total mass
def mass(N,D,L): # the mass of the fins in total
    return N * (density * L * (2*pi * (D**2 / 4)))

# The gaols is to ensure this function is equal to the total heat-loss
def Q(N,D,L): # total heat loss -> just defiing the function
    q_bare = h * pi * (D**2 / 4) * dt
    
    A = h * pi**2 * k * (D**3/4) * dt
    B = math.sinh(m*L) + (h / (m*k))*math.cosh(m*L)
    C = math.cosh(m*L) + (h / (m*k))*math.sinh(m*L)
    
    q_fin = A * (B/C)
    
    return N*(q_bare + q_fin)

# Wrapper function to pass the variables N, D, and L as a single arg
def objective_function(ID):
    N, D, L = ID
    return mass(N, D, L)

# Wrapper function for the constraint
def constraint_equation(ID):
    N, D, L = ID
    total_heat_loss = q
    return Q(N, D, L) - total_heat_loss

# Setting up the constraint dictionary for the minimize function
constraints = {'type': 'eq', 'fun': lambda ID: constraint_equation(ID)}

# Initial guess for the independent variables [N, D, L]
initial_guess = [1, 1, 1]

# Bounds for the variables (e.g., N, D, L should be positive)
bounds = [(1, None), (1E-5, None), (1E-5, None)]

# Perform the optimization
result = minimize(objective_function, initial_guess, constraints=constraints, bounds=bounds, method='SLSQP')

# Print the result
if result.success:
    print(f'Minimum value of mass: {result.fun}')
    print(f'Optimal values of N, D, and L: {result.x}')
else:
    print('Optimization failed:', result.message)