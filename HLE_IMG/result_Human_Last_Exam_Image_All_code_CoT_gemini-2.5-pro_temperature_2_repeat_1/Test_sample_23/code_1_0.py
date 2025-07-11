import math

# Plan:
# 1. Define all given constants.
# 2. Calculate intermediate values like pipe diameter, area, and fluid velocity.
# 3. Calculate the total head loss (h_L) by summing major (friction) and minor (fittings) losses.
# 4. Use the Bernoulli equation between the tank surface and the mold exit to solve for the pump head (h_p).
# 5. Calculate the pump power (work rate) using the pump head.
# 6. Print a detailed breakdown of the calculation.

# 1. Define given constants
rho = 997         # Density of water in kg/m^3
g = 9.81          # Acceleration due to gravity in m/s^2
Q = 2.86e-3       # Volume flow rate in m^3/s
z1 = 2            # Elevation of tank surface relative to the pump in m
z2 = 3            # Elevation of mold exit relative to the pump in m
r = 15.5e-3       # Pipe radius in m
L = 14.9          # Pipe length in m
f = 0.004         # Darcy friction factor
Le_D_bends = 31   # Equivalent length ratio for bends
K_c = 0.4         # Loss coefficient for entrance (shrinkage/contraction)
K_e = 0.8         # Loss coefficient for exit (expansion)
# Note: P1 and P2 are both at ambient pressure, so they cancel out. v1 and v2 are negligible.

# 2. Calculate pipe parameters and flow velocity
D = 2 * r
A = math.pi * r**2
v = Q / A

# Calculate the velocity head, a common term in all loss calculations
velocity_head = v**2 / (2 * g)

# 3. Calculate head losses
# Major loss due to pipe friction
h_L_major = f * (L / D) * velocity_head
# Minor losses from all fittings
# Convert equivalent length ratio for bends to a loss coefficient K
K_bends = f * Le_D_bends
# Total minor loss is sum of entrance, exit, and bends
h_L_minor = (K_c + K_e + K_bends) * velocity_head
# Total head loss
h_L_total = h_L_major + h_L_minor

# 4. Calculate required pump head from the simplified Bernoulli equation: h_p = (z2 - z1) + h_L_total
delta_z = z2 - z1
h_p = delta_z + h_L_total

# 5. Calculate the work of the pump (Power in Watts): W_pump = ρ * Q * g * h_p
work_pump = rho * Q * g * h_p

# 6. Print the detailed breakdown of the calculation
print("This script calculates the work (power) of the pump using the generalized Bernoulli equation.")
print("The overall equation for pump power is W_pump = ρ * Q * g * h_p")
print("\n--- Calculation Steps ---")

print("\nStep 1: Calculate Total Head Loss (h_L)")
print(f"The fluid velocity in the pipe is v = Q / A = {Q} / ({math.pi:.4f}*{r}^2) = {v:.3f} m/s.")
print(f"The corresponding velocity head is v²/2g = {velocity_head:.4f} m.")
print(f"Major Friction Loss (h_L_major) = f*(L/D)*(v²/2g) = {f}*({L}/{D})*{velocity_head:.4f} = {h_L_major:.4f} m.")
print(f"Minor Losses (h_L_minor) = (K_entrance + K_exit + K_bends)*(v²/2g) = ({K_c} + {K_e} + {K_bends:.3f})*{velocity_head:.4f} = {h_L_minor:.4f} m.")
print(f"Total Head Loss (h_L_total) = {h_L_major:.4f} m + {h_L_minor:.4f} m = {h_L_total:.4f} m.")

print("\nStep 2: Calculate Pump Head (h_p)")
print("From the Bernoulli equation, h_p = (z₂ - z₁) + h_L_total")
print(f"h_p = ({z2} m - {z1} m) + {h_L_total:.4f} m = {h_p:.4f} m.")

print("\nStep 3: Calculate the Work (Power) of the Pump")
print(f"The equation is: W_pump = ρ * Q * g * h_p")
print(f"W_pump = {rho} * {Q} * {g} * {h_p:.4f}")
print(f"W_pump = {work_pump:.1f} W")

# The final answer in the required format
final_answer = work_pump
print(f"\n<<<{final_answer:.1f}>>>")