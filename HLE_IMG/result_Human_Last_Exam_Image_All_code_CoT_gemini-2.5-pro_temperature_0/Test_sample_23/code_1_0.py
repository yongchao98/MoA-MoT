import math

# Given values from the problem description and table
rho = 997  # Density of water in kg/m^3
P1 = 105  # Pressure at source in N/m^2
P_out = 1.33e5  # Absolute pressure at manometer (pump outlet) in N/m^2
Q = 2.86e-3  # Volume flow rate in m^3/s
h_tank = 2  # Amount of water in the tank in meters
z_pipe = 3  # Vertical length of the pipe in meters
r = 15.5 / 1000  # Pipe radius in meters
L = 14.9  # Pipe length in meters
f = 0.004  # Friction factor
ef_shrinkage = 0.4  # Minor loss coefficient for entrance (shrinkage)
g = 9.81  # Acceleration due to gravity in m/s^2

# Step 1: Calculate pipe properties and fluid velocity
D = 2 * r  # Pipe diameter
A = math.pi * r**2  # Pipe cross-sectional area
v_pipe = Q / A  # Fluid velocity in the pipe

# Step 2: Calculate head loss in the suction line (from tank to pump inlet)
# This assumes the entire pipe length L is on the suction side.
# Head loss = (Entrance Loss + Pipe Friction Loss)
# h_L = (K_entrance + f * L/D) * v^2 / (2g)
L_over_D = L / D
kinetic_head = v_pipe**2 / (2 * g)
h_L_suction = (ef_shrinkage + f * L_over_D) * kinetic_head

# Step 3: Calculate the pressure at the pump inlet (P_in) using Bernoulli's equation
# P1/(rho*g) + z1 + v1^2/(2g) = Pin/(rho*g) + zin + vin^2/(2g) + h_L_suction
# We set the pump inlet as the datum z_in = 0.
# The elevation of the water surface z1 is h_tank + z_pipe.
z1_minus_zin = h_tank + z_pipe
# v1 is 0 (large tank). vin is v_pipe.
# Rearranging for Pin:
# Pin = rho * g * (P1/(rho*g) + z1_minus_zin - vin^2/(2g) - h_L_suction)
P_in_head = (P1 / (rho * g)) + z1_minus_zin - kinetic_head - h_L_suction
P_in = P_in_head * rho * g

# Step 4: Calculate the work of the pump per unit mass
# W_p = (P_out - P_in) / rho
work_pump = (P_out - P_in) / rho

# Output the final equation with calculated values
print("Calculation Steps:")
print(f"1. Pipe velocity (v) = {v_pipe:.3f} m/s")
print(f"2. Suction head loss (h_L_suction) = {h_L_suction:.3f} m")
print(f"3. Pump inlet pressure (P_in) = {P_in:.2f} N/m^2")
print("\nFinal Equation for Pump Work (W_p):")
print(f"W_p = (P_out - P_in) / ρ")
print(f"W_p = ({P_out:.2f} N/m^2 - {P_in:.2f} N/m^2) / {rho} kg/m^3")
print(f"W_p = {work_pump:.2f} J/kg")

# Final answer in the required format
final_answer = work_pump
# <<<108.10>>>