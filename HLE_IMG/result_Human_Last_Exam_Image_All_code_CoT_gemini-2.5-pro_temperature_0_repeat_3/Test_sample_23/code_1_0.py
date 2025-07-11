import math

# 1. Define given variables and constants
rho = 997  # Density of water in kg/m^3
P_M = 1.33e5  # Absolute pressure at manometer in N/m^2
Q = 2.86e-3  # Volume flow rate in m^3/s
g = 9.81  # Acceleration due to gravity in m/s^2

# From the table
r = 15.5 / 1000  # Pipe radius in m
L = 14.9  # Pipe length in m
f = 0.004  # Darcy friction factor
ef_shrinkage = 0.4 # Minor loss coefficient for pipe entrance (sudden contraction)

# From the problem description
h_tank = 2.0  # Water level in the tank in m
h_pipe_vertical = 3.0 # Vertical length of the pipe in m

# 2. State assumptions and define remaining variables for the energy equation
# Assumption: Ambient pressure is 1.0e5 N/m^2, not 105 N/m^2.
P_1 = 1.0e5  # Ambient pressure at tank surface in N/m^2

# Velocity at tank surface (v_1) is negligible as the tank is large.
v_1 = 0

# Set the pump level as the datum (z_M = 0).
# The water surface is z_1 = h_tank + h_pipe_vertical above the pump.
z_M = 0
z_1 = h_tank + h_pipe_vertical

# 3. Calculate intermediate values
# Pipe diameter and area
D = 2 * r
A = math.pi * r**2

# Velocity in the pipe at the manometer (v_M)
v_M = Q / A

# 4. Calculate the terms of the energy equation (per unit mass, J/kg)

# Pressure term
pressure_term = (P_M - P_1) / rho

# Kinetic energy term
kinetic_term = (v_M**2 - v_1**2) / 2

# Potential energy term
potential_term = g * (z_M - z_1)

# Frictional loss term (W_loss = g * h_L)
# h_L includes major loss (friction) and minor loss (entrance)
# h_L = (f * (L/D) + K_entrance) * v_M^2 / (2*g)
L_over_D = L / D
head_loss = (f * L_over_D + ef_shrinkage) * (v_M**2 / (2 * g))
loss_term = g * head_loss

# 5. Calculate the work of the pump
W_pump = pressure_term + kinetic_term - potential_term + loss_term

# 6. Output the results
print("Calculation Steps:")
print(f"  Pipe Diameter (D): {D:.4f} m")
print(f"  Pipe Area (A): {A:.6f} m^2")
print(f"  Water Velocity (v_M): {v_M:.3f} m/s")
print(f"  Elevation of water surface (z_1): {z_1:.1f} m")
print("-" * 20)
print("Energy Equation Terms (J/kg):")
print(f"  Pressure Term ((P_M - P_1)/ρ): {pressure_term:.3f}")
print(f"  Kinetic Term (v_M²/2): {kinetic_term:.3f}")
print(f"  Potential Term (g(z_M - z_1)): {potential_term:.3f}")
print(f"  Friction Loss Term (W_loss): {loss_term:.3f}")
print("-" * 20)

print("Final Equation:")
# Using potential_term = -g*z_1 since z_M=0
print(f"W_pump = ({P_M:.0f} - {P_1:.0f}) / {rho} + {v_M:.3f}² / 2 + {g} * (0 - {z_1:.1f}) + {loss_term:.3f}")
print(f"W_pump = {pressure_term:.3f} + {kinetic_term:.3f} + ({potential_term:.3f}) + {loss_term:.3f}")
print(f"Work of the pump (W_pump) = {W_pump:.2f} J/kg")

# Final answer in the required format
final_answer = round(W_pump, 2)
print(f"\n<<< {final_answer} >>>")