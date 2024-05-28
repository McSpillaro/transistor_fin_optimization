#%% IMPORTS
import numpy as np
from scipy.optimize import minimize
import json
import logging
from logging_config import setup_logging

#%% LOGGING SETUP
log_file = setup_logging()
logging.info("Optimzation process started.")

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
pi = np.pi # just to shorten it so we don't have to type out the 'np' part
conversion = 1E-3 # conversion from milimeters to meters
H = transistor_details['height']*conversion # (m) height of the transistor
W = transistor_details['width']*conversion # (m) width of the transistor

# SYSTEM PARAMETERS
density = fin_constants['density'] # density of the fin in kg/m3
overdesign_constraint = 0.2 # percent overdesign taken into account
q = system_details['total_heat_loss']*overdesign_constraint # (W) heat loss of the system 
Tw = system_details['surface_temp'] # (Celsius) top surface temperature that must be maintained
T_inf = system_details['environment_temp'] # (Celsius) surroudning environment temperature
dt = Tw - T_inf # (Celsius) temperature difference between surface and environment

# HEAT TRANSFER COEFFICIENT
h = system_details['film_heat_transfer_coefficient'] # (W/m^2*C) film heat transfer coefficient

# THERMCAL CONDUCTIVITY
# Interpolating for the proepr k-value of aluminum at 85 celsius
slope = (fin_constants['thermal_conductivity@100']-fin_constants['thermal_conductivity@0']) / (fin_constants['thermal_temp100'] - fin_constants['thermal_temp0'])
k = slope*(Tw-fin_constants['thermal_temp0']) + fin_constants['thermal_conductivity@0'] # (W/m*C) FINAL THERMAL CONDUCTIVY AT 85 CELSIUS

#%% MIN-MAX CALCULATION
# The goal is to MINIMIZE the total mass
def mass(N, D, L): # total mass of the fin
    return N * (density * L * (2 * np.pi * (D**2 / 4)))

def Q_bare(N, D): # defining the bare surface heat loss function
    A_b = transistor_details['height']*transistor_details['width'] # area of the bare surface
    A_cf = pi*(D**2 / 4) # cross-sectional area of 1 fin taking up the space of the transistor surface
    A_s = A_b - N*A_cf # total surface area of the transistor exposed
    return h * A_s  * dt # heat lost by surface as a func. of fin dimensions

def Q_fin(D, L): # defining the pin fin heat loss function
    m = np.sqrt((4 * h) / (D * k))  # recalculate m based on current D and L
    mL = m * L # defining mL for simplification
    
    # threshold for correcting overflow error - approximation of hyperbolic trig. func.
    if mL > 700:
        sinh_mL = 0.5 * np.exp(mL)
        cosh_mL = 0.5 * np.exp(mL)
    elif mL < -700:
        sinh_mL = -0.5 * np.exp(-mL)
        cosh_mL = 0.5 * np.exp(-mL)
    else:
        sinh_mL = np.sinh(mL)
        cosh_mL = np.cosh(mL)
    
    P = D * pi # perimeter of the pin fin (cylindrical)
    A_f = P * L # surface area of the pin fin (cylindrical)
    # spreading out the equation for heat lost by fin for coding simplification
    numerator = sinh_mL + (h / (m * k)) * cosh_mL
    denominator = cosh_mL + (h / (m * k)) * sinh_mL
    coefficient = h * P * k * A_f * dt
    
    return coefficient * (numerator / denominator) # heat lost by a single fin    

# this is the TOTAL heat loss as a func. of pin fin dimensions : CONSTRAINT
def Q_total(N, D, L): # the total heat lost from the system: both bare and fin
    q_bare = Q_bare(N, D) # heat lost by transistor surface as a func. of fin dimensions
    q_fin = Q_fin(D, L) # heat lost by a single pin fin
    
    return q_bare + (N * q_fin)

# fin effectivness which is what we want to MAXIMIZE
def fin_effectiveness(N, D, L):
    m = np.sqrt((4 * h) / (D * k))  # recalculate m based on current D and L
    mL = m * L # defning for simplification
    
    # threshold for correcting overflow error - approximation of hyperbolic trig. func.
    if mL > 700:
        tanh_mL = 1
    elif mL < -700:
        tanh_mL = -1
    else:
        tanh_mL = np.tanh(mL)
    
    return tanh_mL / m
# upper-bound for fin effectivness defining 100% effectiveness constraint
def fin_effectiveness_upper(ID):
    N, D, L = ID
    effectiveness_value = fin_effectiveness(N, D, L)
    logging.info(f'Fin Effectiveness Upper Constraint - N: {N}, D: {D}, L: {L}, Value: {effectiveness_value}')
    return 2 - effectiveness_value

# lower-bound for fin effectiveness defining 0% effectivness constraint
def fin_effectiveness_lower(ID):
    N, D, L = ID
    effectiveness_value = fin_effectiveness(N, D, L)
    logging.info(f'Fin Effectiveness Lower Constraint - N: {N}, D: {D}, L: {L}, Value: {effectiveness_value}')
    return effectiveness_value - 1

def objective_function(ID): # what is to be MIN-MAX
    N, D, L = ID
    objective_value = mass(N, D, L) - fin_effectiveness(N, D, L) # a max fin_effectiveness will minimize mass
    
    # logging the objective value results
    logging.info(f'Objective Function - N: {N}, D: {D}, L: {L}, Value: {objective_value}')
    return objective_value

def constraint_equation(ID): # the CONSTRAINT
    N, D, L = ID
    constraint_value = Q_total(N, D, L) - q # must always = 0 to ensure Q_total always = q

    # logging the constraint value results
    logging.info(f'Constraint Equation - N: {N}, D: {D}, L: {L}, Value: {constraint_value}')
    return constraint_value

# constraints for the functions to maintain reasonable results
constraints = [
    {'type': 'eq', 'fun': lambda ID: constraint_equation(ID)},
    {'type': 'ineq', 'fun': lambda ID: fin_effectiveness_upper(ID)},
    {'type': 'ineq', 'fun': lambda ID: fin_effectiveness_lower(ID)}
]

initial_guess = [1, 0.01, 0.01]

bounds = [(0, None), (1E-5, None), (1E-5, None)]

# Initializing the vars for storing the best results based on optimization
best_result = None
best_objective_value = float('inf')

# Loop over a range of integer values for N
for N in range(1, 101):  # Assuming a reasonable upper limit for N
    logging.info(f'Trying N = {N}')
    
    result = minimize(objective_function, [N, initial_guess[1], initial_guess[2]], 
                      constraints=constraints, bounds=bounds, method='SLSQP')
    
    if result.success and result.fun < best_objective_value:
        best_result = result
        best_objective_value = result.fun

if best_result:
    logging.info(f'Minimum value of mass - R: {best_result.fun}')
    logging.info(f'Optimal values of N, D, and L: {best_result.x}')
    print(f'Minimum value of mass - R: {best_result.fun}')
    print(f'Optimal values of N, D, and L: {best_result.x}')
else:
    logging.error('Optimization failed: No feasible solution found')
    print('Optimization failed: No feasible solution found')