import math

# Given values from the problem description and table
rho = 997  # Density of water in kg/m^3
g = 9.81  # Acceleration due to gravity in m/s^2
P1 = 105  # Ambient pressure at the tank surface in N/m^2
P_M = 1.33e5  # Absolute pressure at the manometer in N/m^2
Q = 2.86e-3  # Volume flow rate in m^3/s
z1 = 2  # Elevation of water in the tank in m
z2 = 3  # Elevation of the manometer/mold inlet in m (same as z_M)
r = 15.5 / 1000  # Pipe radius in m
L = 14.9  # Pipe length in m
f = 0.004  # Friction factor
K_c = 0.4  # Minor loss coefficient for contraction (e_f shrinkage)

# --- Step 1: Calculate pipe parameters and fluid velocity ---
D = 2 * r  # Pipe diameter in m
A = math.pi * r**2  # Pipe cross-sectional area in m^2
v_M = Q / A  # Velocity of water in the pipe (at the manometer) in m/s

# --- Step 2: Calculate head terms ---
# Velocity at tank surface is negligible (v1 = 0)
v1 = 0

# Velocity head at the manometer
velocity_head_M = v_M**2 / (2 * g)

# --- Step 3: Calculate head losses between tank and manometer ---
# Major loss due to pipe friction
L_over_D = L / D
major_loss = f * L_over_D * velocity_head_M

# Minor loss due to sharp entrance from tank to pipe
minor_loss = K_c * velocity_head_M

# Total head loss
h_losses = major_loss + minor_loss

# --- Step 4: Apply the extended Bernoulli's equation to find the pump head (h_pump) ---
# h_pump = (P_M - P1)/(rho*g) + (v_M^2 - v1^2)/(2*g) + (z2 - z1) + h_losses
pressure_head_diff = (P_M - P1) / (rho * g)
velocity_head_diff = (v_M**2 - v1**2) / (2 * g)
elevation_head_diff = z2 - z1

h_pump = pressure_head_diff + velocity_head_diff + elevation_head_diff + h_losses

# --- Step 5: Calculate the work of the pump per unit mass ---
w_pump = h_pump * g

# --- Step 6: Print the final equation with all values ---
print("Calculation of the work of the pump (w_pump):")
print("\nBased on the energy equation: w_pump = g * h_pump")
print("where h_pump = (P_M - P₁) / (ρ*g) + (v_M² - v₁²) / (2*g) + (z_M - z₁) + h_losses")
print("h_losses = h_major + h_minor = (f * (L/D) + K_c) * (v_M² / (2*g))")
print("\nSubstituting the values into the equation for w_pump:")
# Equation for work per unit mass w_pump = (P_M-P1)/rho + g*(z2-z1) + (v_M^2-v1^2)/2 + g*h_losses
# It's clearer to print the h_pump calculation first.

print(f"h_pump = (({P_M:.0f} - {P1}) / ({rho} * {g})) + (({v_M:.3f}² - {v1}²) / (2 * {g})) + ({z2} - {z1}) + ({f} * ({L}/{D:.3f}) + {K_c}) * (({v_M:.3f}²) / (2 * {g}))")
print(f"h_pump = ({pressure_head_diff:.3f}) + ({velocity_head_diff:.3f}) + ({elevation_head_diff}) + ({h_losses:.3f})")
print(f"h_pump = {h_pump:.3f} m")

print(f"\nw_pump = h_pump * g = {h_pump:.3f} * {g}")
print(f"The calculated work of the pump is: {w_pump:.2f} J/kg")

# Final answer in specified format
print(f"\n<<<{w_pump:.2f}>>>")