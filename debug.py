# main_script.py

import numpy as np
from scipy.optimize import minimize
import json
import logging
from logging_config import setup_logging

# Set up logging
log_file = setup_logging()

logging.info("Optimization process started.")

# FILE HANDLING
with open('data.json') as f:
    data = json.load(f)

transistor_details_RAW = data['transistor_details']
fin_constants_RAW = data['fin_constants']
system_details_RAW = data['system_details']

transistor_details = {}
fin_constants = {}
system_details = {}

for detail in transistor_details_RAW:
    transistor_details.update(detail)
for detail in fin_constants_RAW:
    fin_constants.update(detail)
for detail in system_details_RAW:
    system_details.update(detail)

pi = np.pi
conversion = 1E-3
H = transistor_details['height'] * conversion
W = transistor_details['width'] * conversion

density = fin_constants['density']
overdesign_constraint = 0.2
q = system_details['total_heat_loss'] * overdesign_constraint
Tw = system_details['surface_temp']
T_inf = system_details['environment_temp']

h = system_details['film_heat_transfer_coefficient']

slope = (fin_constants['thermal_conductivity@100'] - fin_constants['thermal_conductivity@0']) / (fin_constants['thermal_temp100'] - fin_constants['thermal_temp0'])
k = slope * (Tw - fin_constants['thermal_temp0']) + fin_constants['thermal_conductivity@0']
dt = Tw - T_inf

def mass(N, D, L):
    return N * (density * L * (2 * np.pi * (D**2 / 4)))

def Q_bare(N, D):
    A_b = transistor_details['height'] * transistor_details['width']
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
    return tanh_mL / m

def objective_function(ID):
    N, D, L = ID
    objective_value = mass(N, D, L) - fin_effectiveness(N, D, L)
    logging.info(f'Objective Function - N: {N}, D: {D}, L: {L}, Value: {objective_value}')
    return objective_value

def constraint_equation(ID):
    N, D, L = ID
    constraint_value = Q_total(N, D, L) - q
    logging.info(f'Constraint Equation - N: {N}, D: {D}, L: {L}, Value: {constraint_value}')
    return constraint_value

# constraints for the functions to maintain reasonable results
constraints = [
    {'type': 'eq', 'fun': lambda ID: constraint_equation(ID)}
    # {'type': 'ineq', 'fun': lambda ID: fin_effectiveness_upper(ID)},
    # {'type': 'ineq', 'fun': lambda ID: fin_effectiveness_lower(ID)}
]

# Initial guesses for [N, D, L]
initial_guess = [1, 0.01, 0.01]
# Boundary conditions for [N, D, L]
bounds = [(0, None), (1E-5, None), (1E-5, None)]

# Initializing the vars for storing the best results based on optimization
best_result = None
best_objective_value = float('inf')

# Loop over a range of integer values for N
for N in range(1, 101):  # Assuming a reasonable upper limit for N
    logging.info(f'Trying N = {N}')
    # THE ACTUAL FUNCTION FOR MINIMIZING ALL FUNCTIONS DEFINED ABOVE
    result = minimize(objective_function, [N, initial_guess[1], initial_guess[2]], 
                      constraints=constraints, bounds=bounds, method='SLSQP')
    # Checks for the best valid result
    if result.success and result.fun < best_objective_value:
        best_result = result
        best_objective_value = result.fun

# Logs results and prints the best result to the terminal
if best_result:
    logging.info(f'Minimum value of mass - R: {best_result.fun}')
    logging.info(f'Optimal values of N, D, and L: {best_result.x}')
    print(f'Minimum value of mass - R: {best_result.fun}')
    print(f'Optimal values of N, D, and L: {best_result.x}')
else:
    logging.error('Optimization failed: No feasible solution found')
    print('Optimization failed: No feasible solution found')