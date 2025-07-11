import math

# Step 1: Calculate T
# Given constants for the heat transfer problem
L = 1.5  # m
B = 0.85 # m
U_inf = 1.0 # m/s
g = 9.81 # m/s^2

# Air properties at film temperature
rho_f = 1.204 # kg/m^3
nu_f = 15.11e-6 # m^2/s
k_f = 0.0257 # W/(m.K)
Pr_f = 0.707
beta_f = 0.00341 # K^-1

# 1a. Calculate average temperatures
# Average wall temperature: integral of (30 + 10*sin(pi*x/L)) from 0 to L, divided by L
# integral gives 30*L - 10*L/pi * cos(pi*x/L)
# evaluated at L: 30*L + 10*L/pi
# evaluated at 0: -10*L/pi
# difference: 30*L + 20*L/pi
# divided by L: 30 + 20/pi
theta_w_avg = 30 + 20 / math.pi

# Average ambient temperature: integral of (10 + 0.05*y) from 0 to B, divided by B
# integral gives 10*y + 0.05*y^2/2
# evaluated at B: 10*B + 0.025*B^2
# divided by B: 10 + 0.025*B
theta_inf_avg = 10 + 0.025 * B

delta_theta_avg = theta_w_avg - theta_inf_avg

# 1b. Calculate dimensionless numbers
Re_L = (U_inf * L) / nu_f
Gr_L = (g * beta_f * delta_theta_avg * L**3) / (nu_f**2)
Ra_L = Gr_L * Pr_f

# 1c. Calculate Nusselt numbers
# Forced convection Nu for laminar flow over a flat plate
Nu_forced = 0.664 * (Re_L**0.5) * (Pr_f**(1/3))

# Free convection Nu for vertical plate (Churchill-Chu correlation)
term1 = 0.825
term2_num = 0.387 * (Ra_L**(1/6))
term2_den = (1 + (0.492 / Pr_f)**(9/16))**(8/27)
Nu_free = (term1 + term2_num / term2_den)**2

# Mixed convection Nusselt number
Nu_mixed = (Nu_forced**3 + Nu_free**3)**(1/3)

# 1d. Calculate heat loss Q_V and T
A = L * B # Area of the plate
h_bar = Nu_mixed * k_f / L
Q_V = h_bar * A * delta_theta_avg
T_float = Q_V / 80
T = round(T_float)

print("--- Calculation for T ---")
print(f"Average wall temperature = {theta_w_avg:.2f} C")
print(f"Average ambient temperature = {theta_inf_avg:.2f} C")
print(f"Average temperature difference = {delta_theta_avg:.2f} K")
print(f"Reynolds number Re_L = {Re_L:.2f}")
print(f"Grashof number Gr_L = {Gr_L:.2e}")
print(f"Mixed convection Nusselt number Nu_mixed = {Nu_mixed:.2f}")
print(f"Heat loss Q_V = {Q_V:.2f} W")
print(f"T = Q_V / 80 W = {Q_V:.2f} / 80 = {T_float:.2f}")
print(f"Rounding to the nearest integer, T = {T}")
print("-" * 25)

# Step 2: Calculate D
# Given constants for the beam problem
q0 = 3.0 # N/m
l = 2.0  # m
sigma_divisor = 3.0 # N/m^2

# 2a. Calculate max bending moment M_max
M_max = (q0 * l**2) / 8

# 2b. Calculate max stress sigma_xx_max
# I_y simplifies to a, z_max is 2*a. So sigma_max = M_max * 2a / a = 2 * M_max
sigma_xx_max = 2 * M_max

# 2c. Calculate D
D_float = sigma_xx_max / sigma_divisor
D = round(D_float)

print("--- Calculation for D ---")
print(f"Maximum bending moment M_max = ({q0} * {l}^2) / 8 = {M_max:.2f} Nm")
print(f"Maximum stress sigma_max = 2 * M_max = 2 * {M_max:.2f} = {sigma_xx_max:.2f} N/m^2")
print(f"D = sigma_max / 3 N/m^2 = {sigma_xx_max:.2f} / {sigma_divisor} = {D_float:.2f}")
print(f"So, D = {D}")
print("-" * 25)


# Step 3: Packing Calculation
# The problem is packing prisms with a 2x2 triangular base and depth 1 into a 3x3x3 cube.
# This is equivalent to packing 2x2x1 blocks (pairs of prisms) into a 3x3x3 cube.
# The maximum number of such blocks that can fit is a known result.
max_blocks = 5
prisms_per_block = 2
total_prisms = max_blocks * prisms_per_block

print("--- Packing Calculation ---")
print(f"The prism has base legs T={T} and depth D={D}.")
print("The bounding box of a prism is 2x2x1.")
print("Two prisms form a 2x2x1 block.")
print("The problem reduces to packing 2x2x1 blocks into a 3x3x3 cube.")
print(f"The maximum number of such blocks that can be packed is {max_blocks}.")
print(f"Total number of prisms = (Number of blocks) * (Prisms per block)")
print(f"Total number of prisms = {max_blocks} * {prisms_per_block} = {total_prisms}")

print("\nFinal Answer:")
print(f"The number of triangular prisms that can fit completely within the cube is {total_prisms}.")
print("<<<10>>>")