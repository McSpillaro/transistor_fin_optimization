import math as m
import numpy as np

# DIMENSIONS OF THE TRANSISTOR SURFACE
conversion = 1E-3 # conversion from milimeters to meters
H = 50*conversion # (m) height of the transistor
W = 100*conversion # (m) width of the transistor

# SYSTEM PARAMETERS
q = 50 # (W) heat loss of the system 
T0 = 85 # (Celsius) top surface temperature that must be maintained
T1 = 30 # (Celsius) surroudning environment temperature

# HEAT TRANSFER COEFFICIENT
h = 15 # (W/m^2*C) film heat transfer coefficient

# THERMCAL CONDUCTIVITY
# Interpolating for the proepr k-value of aluminum at 85 celsius
kT0 = 0 # (Celsius) temperature used for k0 in Appendix A.2
kT1 = 100 # (Celsius) temperature used for k1 in Appendix A.2
k0 = 202 # (W/m*C) thermal conducitivy @ 0 celsius
k1 = 206 # (W/m*C) thermal conductivity @ 100 celsius
mk = (k1-k0) / (kT1-kT0) # slope of the interpolation function
k = mk*(T0-kT0) + k0 # (W/m*C) FINAL THERMAL CONDUCTIVY AT 85 CELSIUS

# m-PARAMETER
L = 1 # (m) initiailzing the length of the pin fin
D = 1 # (m) initializing the diamter of the pin fin
r = D/2 # radius of the pin fin
P = D*m.pi # perimeter (AKA: circumference) of pin fin
Ac = m.pi*r**2 # cross-sectional area of the pin fin
numerator = h*P # numerator specific to the m-parameter
denominator = k*Ac
m = m.sqrt(numerator/denominator)