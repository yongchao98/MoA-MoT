import math

# Given values from the problem description and table
rho = 997  # Density of water in kg/m^3
P1 = 1.05e5  # Ambient pressure at tank surface in N/m^2
Pm = 1.33e5  # Absolute pressure at manometer in N/m^2
Q = 2.86e-3  # Volume flow rate in m^3/s
z1 = 2  # Height of water in the tank in m
zm = 0  # Elevation of the manometer (datum)
g = 9.81  # Acceleration due to gravity in m/s^2

# Pipe and loss properties
r = 15.5 / 1000  # Pipe radius in m
L = 14.9  # Pipe length in m
f = 0.004  # Darcy friction factor
Kent = 0.4  # Minor loss coefficient for entrance (shrinkage)
Le_D_bend = 31 # Equivalent length ratio for bend

# --- Step-by-step calculation ---

# 1. Calculate pipe diameter and area
D = 2 * r
A = math.pi * r**2

# 2. Calculate flow velocity
v = Q / A

# 3. Calculate velocity head term (v^2 / 2)
v_sq_over_2 = v**2 / 2

# 4. Calculate head loss terms multiplied by g (to get work per mass)
# Work loss from entrance
work_loss_ent = Kent * v_sq_over_2
# Work loss from pipe friction
work_loss_friction = f * (L / D) * v_sq_over_2
# Work loss from bend
work_loss_bend = f * Le_D_bend * v_sq_over_2

# Total work loss
total_work_loss = work_loss_ent + work_loss_friction + work_loss_bend

# 5. Calculate terms in the work equation for the pump
# Work from pressure change
work_pressure = (Pm - P1) / rho
# Work from kinetic energy change (v1 is 0)
work_kinetic = v_sq_over_2
# Work from potential energy change
work_potential = g * (zm - z1)

# 6. Calculate total work of the pump per unit mass
work_pump = work_pressure + work_kinetic + work_potential + total_work_loss

# --- Print the results and the final equation ---
print("Calculation of the Work of the Pump (per unit mass)")
print("-" * 50)
print(f"The work of the pump (w_p) is calculated using the energy balance equation:")
print("w_p = (P_m - P_1)/ρ + (v^2)/2 + g*(z_m - z_1) + w_loss\n")
print("Where w_loss is the work to overcome friction losses:\n")
print(f"w_loss = (K_ent + f*(L/D) + f*(L_e/D)_bend) * v^2/2\n")

print("Substituting the values into the equation:")
print(f"w_p = ({Pm:.0f} - {P1:.0f})/{rho} + ({v:.2f}^2)/2 + {g}*({zm} - {z1}) + ({Kent} + {f}*{L}/{D:.3f} + {f}*{Le_D_bend}) * ({v:.2f}^2)/2")

# Print the value of each term in the equation
print("\nCalculating each term:")
print(f"Pressure term = {work_pressure:.2f} J/kg")
print(f"Kinetic energy term = {work_kinetic:.2f} J/kg")
print(f"Potential energy term = {work_potential:.2f} J/kg")
print(f"Friction loss term = {total_work_loss:.2f} J/kg")

# Print the final result
print("\nFinal Calculation:")
print(f"w_p = {work_pressure:.2f} + {work_kinetic:.2f} + ({work_potential:.2f}) + {total_work_loss:.2f}")
print(f"w_p = {work_pump:.2f} J/kg")
print(f"\n<<<The work of the pump is {work_pump:.2f} J/kg.>>>")
print(f'<<<{work_pump:.1f}>>>')