import math

# 1. Define knowns from the problem statement and table
rho = 997  # Density of water (kg/m^3)
P_ambient = 105  # Ambient pressure (N/m^2)
Q = 2.86e-3  # Volume flow rate (m^3/s)
z1 = 2  # Height of water in the tank (m), relative to the pump
z2 = 3  # Vertical length of the pipe from the pump (m)
r = 15.5 / 1000  # Pipe radius (m)
L_div_D = 31  # L/D ratio from the table
f = 0.004  # Friction factor
K_L_entrance = 0.4  # Loss coefficient for shrinkage (entrance)
K_L_exit = 0.8  # Loss coefficient for expansion (exit)
g = 9.81  # Acceleration due to gravity (m/s^2)

# Assumptions based on the problem statement:
# P1 = P2 = P_ambient, so the pressure head term (P2-P1)/rho*g is 0.
# v1 (tank surface velocity) is negligible (v1=0).
# v2 (exit velocity) kinetic energy is negligible (v2=0).
# The Bernoulli equation simplifies to: z1 + H_p = z2 + h_L
# So, H_p = (z2 - z1) + h_L

# 2. Calculate intermediate values
D = 2 * r  # Pipe diameter (m)
A = math.pi * r**2  # Pipe cross-sectional area (m^2)
v = Q / A  # Water velocity in the pipe (m/s)
v_head = v**2 / (2 * g) # Velocity head (m)

# 3. Calculate Head Losses (h_L)
# Major head loss
h_L_major = f * L_div_D * v_head
# Minor head loss
h_L_minor = (K_L_entrance + K_L_exit) * v_head
# Total head loss
h_L_total = h_L_major + h_L_minor

# 4. Calculate Pump Head (H_p)
delta_z = z2 - z1
H_p = delta_z + h_L_total

# 5. Calculate Work of the Pump (w_pump)
w_pump = H_p * g

# 6. Print the results and the final equation
print("--- Calculation Steps ---")
print(f"1. Pipe velocity (v) = Q / A = {Q:.4f} / ({math.pi:.3f} * {r:.4f}^2) = {v:.2f} m/s")
print(f"2. Velocity head (v^2 / 2g) = {v:.2f}^2 / (2 * {g}) = {v_head:.3f} m")
print("\n3. Head Losses:")
print(f"   - Major Loss (h_L_major) = f * (L/D) * (v^2 / 2g) = {f} * {L_div_D} * {v_head:.3f} = {h_L_major:.3f} m")
print(f"   - Minor Loss (h_L_minor) = (K_entrance + K_exit) * (v^2 / 2g) = ({K_L_entrance} + {K_L_exit}) * {v_head:.3f} = {h_L_minor:.3f} m")
print(f"   - Total Head Loss (h_L_total) = {h_L_major:.3f} + {h_L_minor:.3f} = {h_L_total:.3f} m")
print(f"\n4. Change in Elevation Head (z2 - z1) = {z2} - {z1} = {delta_z:.3f} m")
print(f"\n5. Total Pump Head (H_p) = (z2 - z1) + h_L_total = {delta_z:.3f} + {h_L_total:.3f} = {H_p:.3f} m")

print("\n--- Final Equation and Answer ---")
print("The work of the pump per unit mass (w_pump) is calculated as:")
print("w_pump = g * H_p = g * ((z2 - z1) + h_L_total)")
print(f"w_pump = {g} * (({z2} - {z1}) + {h_L_total:.3f})")
print(f"w_pump = {g} * ({H_p:.3f})")
print(f"w_pump = {w_pump:.2f} J/kg")

<<<19.31>>>