import math

# 1. Define all given parameters
rho = 997  # Density of water in kg/m^3
g = 9.81  # Acceleration due to gravity in m/s^2
Q = 2.86e-3  # Volume flow rate in m^3/s

# Heights (assuming pump is at z=0)
z1 = 2.0  # Height of water in the tank in meters
z2 = 3.0  # Vertical length of the discharge pipe in meters

# Pipe and loss parameters from the table
r_mm = 15.5  # Pipe radius in mm
length = 14.9  # Pipe length in m
L_over_D = 31  # L/D ratio for friction calculation
f = 0.004  # Friction factor
ef_shrinkage = 0.4  # Minor loss coefficient for entrance (shrinkage)
ef_expansion = 0.8  # Minor loss coefficient for exit (expansion)

# Note on inconsistent data: The provided r and length give L/D = 14.9 / (2*0.0155) = 480.6.
# This contradicts the table value L/D = 31. We will proceed using the explicit L/D value from the table as intended for the friction loss calculation.
# The manometer pressure is not needed for this solution method, which relies on an overall energy balance.

# 2. Calculate intermediate values
r = r_mm / 1000.0  # Convert radius to meters
D = 2 * r  # Pipe diameter in meters
A = math.pi * r**2  # Pipe cross-sectional area in m^2
v_pipe = Q / A  # Velocity of water in the pipe in m/s
v_sq_over_2 = (v_pipe**2) / 2 # Velocity head term in J/kg (m^2/s^2)

# 3. Calculate the work lost to friction (W_loss)
# W_loss = (major loss term + minor loss term) * velocity_head
K_total = ef_shrinkage + ef_expansion
W_loss = (f * L_over_D + K_total) * v_sq_over_2

# 4. Calculate the change in potential energy (W_potential)
W_potential = g * (z2 - z1)

# 5. Calculate the total work of the pump (W_pump)
# Based on the simplified Bernoulli's equation: W_pump = W_potential + W_loss
# Assumptions: P1=P2, v1=0, v2=0
W_pump = W_potential + W_loss

# 6. Print the results step-by-step
print("--- Calculation Steps ---")
print(f"Pipe velocity (v): {v_pipe:.3f} m/s")
print(f"Velocity head (v^2 / 2): {v_sq_over_2:.3f} J/kg")
print("\nWork of the pump is the sum of potential energy change and friction losses.")
print("Equation: W_pump = g * (z2 - z1) + (f * (L/D) + K_total) * (v^2 / 2)")
print("\n--- Substituting Values ---")
print(f"W_pump = {g} * ({z2} - {z1}) + ({f} * {L_over_D} + ({ef_shrinkage} + {ef_expansion})) * ({v_pipe:.3f}^2 / 2)")
print(f"W_pump = {g} * ({z2-z1}) + ({f * L_over_D:.3f} + {K_total}) * {v_sq_over_2:.3f}")
print(f"W_pump = {W_potential:.3f} J/kg (for potential energy) + {W_loss:.3f} J/kg (for losses)")
print(f"W_pump = {W_pump:.3f} J/kg")

print("\n--- Final Answer ---")
print(f"The calculated work of the pump is {W_pump:.2f} J/kg.")
print(f"The final equation with all numbers is:")
print(f"Work of Pump = ({g} * ({z2} - {z1})) + (({f} * {L_over_D}) + {ef_shrinkage} + {ef_expansion}) * (({v_pipe:.3f}**2) / 2)")
print(f"Work of Pump = {W_potential:.2f} + ({f * L_over_D:.3f} + {K_total}) * {v_sq_over_2:.2f}")
print(f"Work of Pump = {W_potential:.2f} + {W_loss:.2f} = {W_pump:.2f} J/kg")
<<<19.32>>>