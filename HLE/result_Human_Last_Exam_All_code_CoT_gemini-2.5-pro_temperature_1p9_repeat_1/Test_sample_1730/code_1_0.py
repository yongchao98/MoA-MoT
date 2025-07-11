import math

# Part 1: Calculate T
# --------------------
# Given constants for heat transfer
L = 1.5  # m, length of the collector
B = 0.85  # m, width of the collector
U_inf = 1.0  # m/s, wind speed
g = 9.81  # m/s^2

# Air properties at film temperature
k_f = 0.0257  # W/(m.K), thermal conductivity
nu_f = 15.11e-6  # m^2/s, kinematic viscosity
Pr_f = 0.707  # Prandtl number
beta_f = 0.00341  # K^-1, thermal expansion coefficient

# 1a. Calculate average temperatures
# Average surface temperature: integral of (30 + 10*sin(pi*x/L)) from 0 to L, divided by L
# integral of sin(pi*x/L) dx = [-L/pi * cos(pi*x/L)] from 0 to L = 2L/pi
theta_w_avg = 30.0 + 10.0 * (2.0 / math.pi)

# Average ambient temperature: integral of (10 + 0.05y) from 0 to B, divided by B
theta_inf_avg = 10.0 + 0.05 * B / 2.0

delta_T_avg = theta_w_avg - theta_inf_avg

# 1b. Calculate Reynolds and Grashof numbers
Re_L = (U_inf * L) / nu_f
Gr_L = (g * beta_f * delta_T_avg * (L**3)) / (nu_f**2)
Ra_L = Gr_L * Pr_f

# 1c. Calculate Nusselt numbers for mixed convection
# Forced convection (laminar flow over flat plate)
Nu_forced = 0.664 * (Re_L**0.5) * (Pr_f**(1/3))

# Natural convection (vertical plate, Churchill-Chu correlation)
term1 = 0.387 * (Ra_L**(1/6))
term2 = (1 + (0.492 / Pr_f)**(9/16))**(8/27)
Nu_natural = (0.825 + term1 / term2)**2

# Mixed convection (transverse flow)
Nu_mixed = (Nu_forced**3 + Nu_natural**3)**(1/3)

# 1d. Calculate total heat loss and T
h_avg = (Nu_mixed * k_f) / L
Area = L * B
Q_V = h_avg * Area * delta_T_avg
T_value = Q_V / 80
T = round(T_value)

# Part 2: Calculate D
# --------------------
# Given constants for beam mechanics
q0 = 3.0  # N/m, distributed load
l_beam = 2.0  # m, beam length

# 2a. Calculate max bending moment for a simply supported beam
M_max = (q0 * l_beam**2) / 8.0

# 2b. Calculate geometrical properties. The problem simplifies to sigma_max = 2 * M_max
# The constant 'a' is not needed for the final calculation of sigma_max due to simplification.
# a = ( (64/3) - (math.pi/4) )**(-1/3)
# I_y = ( (64/3) - (math.pi/4) ) * a**4 = (1/a**3) * a**4 = a
# z_max = 2*a
# sigma_max = (M_max * z_max) / I_y = (M_max * 2*a) / a = 2 * M_max
sigma_max = 2 * M_max

# 2c. Calculate D
D_value = sigma_max / 3.0
D = int(D_value)

# Part 3: Solve Packing Problem
# ------------------------------
# The prism has a right triangular base with legs of length T=2 and depth D=1.
# The cube is 3x3x3.

# 3a. We can stack prisms along the D-dimension.
# Number of layers = Cube side / D
num_layers = 3 / D

# 3b. In one 3x3 layer, we can fit 3 prisms.
# (Two prisms form a 2x2 square in one corner, a third prism fits diagonally opposite)
prisms_per_layer = 3

# 3c. Calculate total number of prisms
total_prisms = prisms_per_layer * num_layers

print(f"The calculated value for T is {T}.")
print(f"The calculated value for D is {D}.")
print("The base of each prism is a right triangle with two sides of length T.")
print(f"The depth of each prism is D={D}.")
print(f"The cube has dimensions 3x3x3.")
print(f"The number of prisms per layer in a 3x3 area is {prisms_per_layer}.")
print(f"The number of available layers in the cube is {int(num_layers)}.")
print(f"The total number of prisms is calculated by: {prisms_per_layer} * {int(num_layers)} = {int(total_prisms)}")
print(f"<<<{int(total_prisms)}>>>")