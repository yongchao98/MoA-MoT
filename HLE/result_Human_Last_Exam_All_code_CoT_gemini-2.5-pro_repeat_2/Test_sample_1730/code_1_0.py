import numpy as np

# --- Part 1: Calculate T ---
print("Part 1: Calculating the prism base side T")

# Given parameters for the solar collector
L = 1.5  # m, Length
B = 0.85  # m, Width
U_inf = 1.0  # m/s, Wind speed
# Air properties at film temperature
k_f = 0.0257  # W/(m.K), Thermal conductivity
nu_f = 15.11e-6  # m^2/s, Kinematic viscosity
Pr_f = 0.707  # Prandtl number
beta_f = 0.00341  # K^-1, Thermal expansion coefficient
g = 9.81  # m/s^2

# Calculate average wall temperature by integrating θ_w(x) over L
# θ_w_avg = (1/L) * ∫[30 + 10*sin(πx/L)]dx from 0 to L = 30 + 20/π
theta_w_avg = 30 + 20 / np.pi
print(f"Average wall temperature: {theta_w_avg:.2f} C")

# Calculate average ambient temperature by integrating θ_∞(y) over B
# θ_∞_avg = (1/B) * ∫[10 + 0.05y]dy from 0 to B = 10 + 0.025*B
theta_inf_avg = 10 + 0.025 * B
print(f"Average ambient temperature: {theta_inf_avg:.2f} C")

# Average temperature difference
delta_theta_avg = theta_w_avg - theta_inf_avg
print(f"Average temperature difference: {delta_theta_avg:.2f} K")

# Calculate dimensionless numbers to determine convection type
Re_L = U_inf * L / nu_f
Gr_L = g * beta_f * delta_theta_avg * L**3 / nu_f**2
print(f"Reynolds number (Re_L): {Re_L:.2e}")
print(f"Grashof number (Gr_L): {Gr_L:.2e}")

# The ratio Gr_L / Re_L^2 is O(1), indicating mixed convection.
# Calculate Nusselt number for forced convection (laminar flow, Re < 5e5)
Nu_L_forced = 0.664 * Re_L**0.5 * Pr_f**(1/3)
print(f"Forced convection Nusselt number: {Nu_L_forced:.1f}")

# Calculate Nusselt number for natural convection (Churchill and Chu correlation)
Ra_L = Gr_L * Pr_f
term1 = 0.825
term2_num = 0.387 * Ra_L**(1/6)
term2_den = (1 + (0.492 / Pr_f)**(9/16))**(8/27)
Nu_L_natural = (term1 + term2_num / term2_den)**2
print(f"Natural convection Nusselt number: {Nu_L_natural:.1f}")

# Combine for mixed convection Nusselt number (with n=3)
Nu_L = (Nu_L_forced**3 + Nu_L_natural**3)**(1/3)
print(f"Combined Nusselt number: {Nu_L:.1f}")

# Calculate average heat transfer coefficient and total heat loss
h_avg = Nu_L * k_f / L
A = L * B
Q_V = h_avg * A * delta_theta_avg
print(f"Total heat loss (Q_V): {Q_V:.2f} W")

# Calculate T and round to the nearest integer
T_float = Q_V / 80.0
T = int(round(T_float))
print(f"Prism base side T = Q_V / 80 = {T_float:.2f}, which is rounded to {T}")
print("-" * 30)

# --- Part 2: Calculate D ---
print("Part 2: Calculating the prism depth D")

# Given parameters for the beam
q0 = 3.0  # N/m, Uniformly distributed load
l = 2.0  # m, Length of the beam

# Calculate maximum bending moment for a simply supported beam
M_max = q0 * l**2 / 8
print(f"Maximum bending moment (M_max): {M_max:.2f} N.m")

# The cross-section properties (I_yy and z_max) are defined in terms of 'a'.
# The specific definition of 'a' implies I_yy = a. The maximum distance z_max = 2a.
# The variable 'a' cancels out in the stress calculation.
# sigma_max = M_max * z_max / I_yy = M_max * (2*a) / a = 2 * M_max
sigma_xx_max = 2 * M_max
print(f"Maximum normal stress (sigma_xx_max): {sigma_xx_max:.2f} N/m^2")

# Calculate D and round to the nearest integer
D_float = sigma_xx_max / 3.0
D = int(round(D_float))
print(f"Prism depth D = sigma_xx_max / 3 = {D_float:.2f}, which is rounded to {D}")
print("-" * 30)

# --- Part 3: Packing Calculation ---
print("Part 3: Calculating the number of prisms in the cube")
cube_dim = 3
prism_base_leg = T
prism_depth = D

print(f"Cube dimensions: {cube_dim} x {cube_dim} x {cube_dim}")
print(f"Prism dimensions: Right triangular base with legs {prism_base_leg}x{prism_base_leg}, depth {prism_depth}")

# The number of layers of prisms that can be stacked is the integer division
# of the cube's dimension by the prism's depth.
num_layers = cube_dim // prism_depth
print(f"Number of layers that can be stacked: floor({cube_dim} / {prism_depth}) = {num_layers}")

# For each layer, we must fit the 2x2 triangular bases into a 3x3 square.
# The most efficient packing allows for 2 triangles per layer.
prisms_per_layer = 2
print(f"Maximum prisms per {cube_dim}x{cube_dim} layer: {prisms_per_layer}")

# The total number of prisms is the product of layers and prisms per layer.
total_prisms = prisms_per_layer * num_layers

print("\nFinal Calculation:")
print("The total number of prisms is (prisms per layer) multiplied by (number of layers).")
print(f"The final equation is: {prisms_per_layer} * {num_layers} = {total_prisms}")

print(f"\n<<< {total_prisms} >>>")