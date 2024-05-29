#%% IMPORTS
import numpy as np
from scipy.optimize import minimize
import json
import logging
from logging_config import setup_logging

#%% LOGGING SETUP
CASE = 3
log_file = setup_logging()
logging.basicConfig(level=logging.INFO, filename=log_file)
logging.info("Optimization process started.")
logging.info(f"CASE: {CASE}")

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

#%% SYSTEM CALCULATIONS FUNCTIONS
def mass(N, D, L): # calculates the total mass of the aluminum fins
    m_value = N * density * L * (pi * (D**2 / 4)) # actual mass value
    logging.debug(f'Calculated mass: {m_value} for N={N}, D={D}, L={L}') # logging
    return m_value

def Q_bare(N, D): # calculates the heat loss of the 
    A_b = H * W # area of the bare surface of the electronic device
    A_cf = pi * (D**2 / 4) # cross-sectional area of the fin in contact with device
    A_void = 1E-4 # area of space between each fin
    A_s = A_b - N * A_cf + A_void # total surface area exposed of the electronic device
    return h * A_s * dt

def Q_fin(D, L): # the heat loss of only one fin
    P = D * pi # perimeter of the fin
    A_f = P * L # surface area of the fin
    
    m = np.sqrt((h * P) / (k * A_f)) # m-value
    mL = m * L # simplification of m * L for ease of typing
    
    n_f = np.tanh(mL) / mL
    return n_f * A_f * h * dt
    
def Q_total(N, D, L): # total heat loss of the fin and device together
    q_bare = Q_bare(N, D)
    q_fin = Q_fin(D, L)
    return q_bare + (N * q_fin)

#%% OPTIMIZATION FUNCTIONS
def fin_effectiveness(N, D, L): # how effective the fin actually is
    P = D * pi # perimeter of the fin
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

def objective_function(ID): # function to optimize mass and fin effectiveness
    N, D, L = ID
    # Add an absolute value to ensure positive mass contribution
    obj_value = abs(mass(N, D, L)) - fin_effectiveness(N, D, L)
    logging.debug(f'Objective function value: {obj_value} for N={N}, D={D}, L={L}')
    return obj_value

def constraint_equation(ID): # function for ensuring calculated heat loss = desired heat loss
    N, D, L = ID
    return Q_total(N, D, L) - q # must be = 0 for constraint to be true

# ensure that mass is always positive
def mass_constraint(ID):
    N, D, L = ID
    m_value = mass(N, D, L) # calculates the mass
    logging.debug(f'Mass constraint value: {m_value} for N={N}, D={D}, L={L}')
    return m_value - 1E-10  # ensure mass is greater than a small positive value

# constraints for the optimization to ensure valid results based on desires
def optimize(nMAX, dMAX, lMAX):
    results = [] # initialize results list
    constraints = [
        {'type': 'eq', 'fun': constraint_equation}, # q_total (calculated) approx. = q (defined)
        {'type': 'ineq', 'fun': mass_constraint}, # mass must be : 0 < mass
    ]

    initial_guess = [1, 0.01, 0.01] # initial guess values for [N, D, L]
    # defined boundaries for [N, D, L]
    bounds = [(1, nMAX), # number of fin boundaries
            (1E-5, dMAX), # boundary for diameter
            (1E-5, lMAX)] # boundary for length

    best_result = None # initializing the best result
    best_objective_value = float('inf') # initializing the best values

    for N in range(1, 101): # loops through tested ranges of number of fins to optimize dimensions
        logging.info(f'Trying N = {N}') # logging
        result = minimize(objective_function, # function to optimize
                        [N, initial_guess[1], initial_guess[2]], # guesses to go into objective func. 
                        constraints=constraints, bounds=bounds,
                        method='SLSQP') # method of optimization
        if result.success and result.fun < best_objective_value: # defines best values
            best_result = result
            best_objective_value = result.fun

    # logging and outputing results
    if best_result:
        logging.info(f'Minimum value of mass: {best_result.fun}')
        logging.info(f'Optimal values of N, D, and L: {best_result.x}')
        results.append(f'Minimum value of mass: {round(best_result.fun*1000, 3)} grams')
        results.append(f'Optimal values of N, D, and L: {best_result.x} -> [# Fins | Diameter (m) | Length (m)]')
    else:
        logging.error('Optimization failed: No feasible solution found')
        print('Optimization failed: No feasible solution found')
        
    return print(f"\n{results[0]}\n{results[1]}\n")

optimize(100, None, 0.0078)