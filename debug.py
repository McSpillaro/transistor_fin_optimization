# file: optimization_script.py

#%% IMPORTS
import numpy as np
from scipy.optimize import minimize
import json
import logging
from logging_config import setup_logging

#%% LOGGING SETUP
log_file = setup_logging()
logging.basicConfig(level=logging.INFO, filename=log_file)
logging.info("Optimization process started.")

#%% FILE HANDLING
# Opening JSON files
with open('data.json') as f:
    data = json.load(f)
    
# Organizing sets of data
transistor_details = {key: val for detail in data['transistor_details'] for key, val in detail.items()}
fin_constants = {key: val for detail in data['fin_constants'] for key, val in detail.items()}
system_details = {key: val for detail in data['system_details'] for key, val in detail.items()}

#%% DEFINING KNOWN SYSTEM VARIABLES
# DIMENSIONS OF THE TRANSISTOR SURFACE
pi = np.pi
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

# THERMAL CONDUCTIVITY
slope = (fin_constants['thermal_conductivity@100'] - fin_constants['thermal_conductivity@0']) / (fin_constants['thermal_temp100'] - fin_constants['thermal_temp0'])
k = slope * (Tw - fin_constants['thermal_temp0']) + fin_constants['thermal_conductivity@0'] # (W/m*C) FINAL THERMAL CONDUCTIVITY AT 85 CELSIUS

#%% MIN-MAX CALCULATION
def mass(N, D, L):
    m_value = N * density * L * (pi * (D**2 / 4))
    logging.debug(f'Calculated mass: {m_value} for N={N}, D={D}, L={L}')
    return m_value

def Q_bare(N, D):
    A_b = H * W
    A_cf = pi * (D**2 / 4)
    A_s = A_b - N * A_cf
    return h * A_s * dt

def Q_fin(D, L):
    m = np.sqrt((4 * h) / (D * k))
    mL = m * L
    if mL > 700:
        sinh_mL = 0.5 * np.exp(mL)
        cosh_mL = 0.5 * np.exp(mL)
    elif mL < -700:
        sinh_mL = -0.5 * np.exp(-mL)
        cosh_mL = 0.5 * np.exp(-mL)
    else:
        sinh_mL = np.sinh(mL)
        cosh_mL = np.cosh(mL)
    P = D * pi
    A_f = P * L
    numerator = sinh_mL + (h / (m * k)) * cosh_mL
    denominator = cosh_mL + (h / (m * k)) * sinh_mL
    coefficient = h * P * k * A_f * dt
    return coefficient * (numerator / denominator)

def Q_total(N, D, L):
    q_bare = Q_bare(N, D)
    q_fin = Q_fin(D, L)
    return q_bare + (N * q_fin)

def fin_effectiveness(N, D, L):
    m = np.sqrt((4 * h) / (D * k))
    mL = m * L
    if mL > 700:
        tanh_mL = 1
    elif mL < -700:
        tanh_mL = -1
    else:
        tanh_mL = np.tanh(mL)
    return tanh_mL / mL

def objective_function(ID):
    N, D, L = ID
    # Add an absolute value to ensure positive mass contribution
    obj_value = abs(mass(N, D, L)) - fin_effectiveness(N, D, L)
    logging.debug(f'Objective function value: {obj_value} for N={N}, D={D}, L={L}')
    return obj_value

def constraint_equation(ID):
    N, D, L = ID
    return Q_total(N, D, L) - q

# Ensure that mass is always positive
def mass_constraint(ID):
    N, D, L = ID
    m_value = mass(N, D, L)
    logging.debug(f'Mass constraint value: {m_value} for N={N}, D={D}, L={L}')
    return m_value - 0.001  # Ensure mass is greater than a small positive value

# New constraint: Fins must be no longer than 1 cm
def length_constraint(ID):
    N, D, L = ID
    return 0.01 - L  # Ensure length is less than or equal to 0.01 meters (1 cm)

constraints = [
    {'type': 'eq', 'fun': constraint_equation},
    {'type': 'ineq', 'fun': mass_constraint},
    {'type': 'ineq', 'fun': length_constraint}
]

initial_guess = [1, 0.01, 0.01]
bounds = [(1, 100), (1E-5, None), (1E-5, 0.01)]

best_result = None
best_objective_value = float('inf')

for N in range(1, 101):
    logging.info(f'Trying N = {N}')
    result = minimize(objective_function,
                      [N, initial_guess[1], initial_guess[2]],
                      constraints=constraints, bounds=bounds,
                      method='SLSQP')
    if result.success and result.fun < best_objective_value:
        best_result = result
        best_objective_value = result.fun

if best_result:
    logging.info(f'Minimum value of mass: {best_result.fun}')
    logging.info(f'Optimal values of N, D, and L: {best_result.x}')
    print('\n')
    print(f'Minimum value of mass: {best_result.fun}')
    print(f'Optimal values of N, D, and L: {best_result.x}')
    print('\n')
else:
    logging.error('Optimization failed: No feasible solution found')
    print('Optimization failed: No feasible solution found')
