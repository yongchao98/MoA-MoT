import math

# 1. Define given values from the problem description and table
rho = 997  # Density of water in kg/m^3
P_ambient = 1.0e5 # Ambient pressure in N/m^2 (Note: this value cancels out)
Q = 2.86e-3  # Volume flow rate in m^3/s
z1 = 2.0  # Elevation of water surface in the tank in meters (relative to pump)
z2 = 3.0  # Elevation of the mold outlet in meters (relative to pump)
r_pipe = 15.5 / 1000  # Pipe radius in meters
L_pipe = 14.9  # Pipe length in meters
f = 0.004  # Darcy friction factor
K_shrinkage = 0.4  # Minor loss coefficient for entrance (shrinkage)
K_expansion = 0.8  # Minor loss coefficient for exit (expansion)
g = 9.81  # Acceleration due to gravity in m/s^2

# 2. Calculate intermediate values
D_pipe = 2 * r_pipe  # Pipe diameter
A_pipe = math.pi * r_pipe**2  # Pipe cross-sectional area
v_pipe = Q / A_pipe  # Fluid velocity in the pipe

# 3. Calculate head losses (h_L)
# The velocity head term v^2 / (2g) is common to all loss calculations
velocity_head = v_pipe**2 / (2 * g)

# Major loss due to pipe friction
h_L_major = f * (L_pipe / D_pipe) * velocity_head

# Minor losses due to fittings (entrance and exit)
K_total_minor = K_shrinkage + K_expansion
h_L_minor = K_total_minor * velocity_head

# Total head loss
h_L_total = h_L_major + h_L_minor

# 4. Calculate the required pump head (h_p)
# From the simplified Bernoulli equation: h_p = (z2 - z1) + h_L_total
delta_z = z2 - z1
h_p = delta_z + h_L_total

# 5. Calculate the work of the pump per unit mass (w_p)
work_pump = g * h_p

# 6. Print the results and the final equation
print("--- Calculation Breakdown ---")
print(f"Pipe velocity (v): {v_pipe:.3f} m/s")
print(f"Velocity head (v^2/2g): {velocity_head:.3f} m")
print(f"Major head loss (friction): {h_L_major:.3f} m")
print(f"Minor head loss (fittings): {h_L_minor:.3f} m")
print(f"Total head loss (h_L): {h_L_total:.3f} m")
print(f"Elevation change (z2 - z1): {delta_z:.3f} m")
print(f"Total pump head required (h_p): {h_p:.3f} m")
print("\n--- Final Calculation for Pump Work ---")
print("The work of the pump (per unit mass) is calculated as: W = g * ( (z2 - z1) + h_L_total )")
print(f"W = {g} * ( ({z2} - {z1}) + {h_L_total:.3f} )")
print(f"W = {g} * ( {delta_z} + {h_L_total:.3f} )")
print(f"W = {g} * ( {h_p:.3f} )")
print(f"W = {work_pump:.2f} J/kg")

<<<32.23>>>