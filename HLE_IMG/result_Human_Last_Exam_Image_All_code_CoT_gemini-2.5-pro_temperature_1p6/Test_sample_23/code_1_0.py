import math

# --- 1. Define Given Variables ---
rho = 997  # Density of water in kg/m^3
g = 9.81  # Acceleration due to gravity in m/s^2
Q = 2.86e-3  # Volume flow rate in m^3/s

# Ambient pressure. The prompt value of 105 N/m^2 is likely a typo for 1.05e5 N/m^2.
# Since P1 = P2, the value cancels out, but we need it for the manometer calculation.
P_ambient = 1.05e5  # Ambient pressure in N/m^2 (Pa)
P_manometer = 1.33e5 # Manometer absolute pressure in N/m^2 (Pa)

# Elevations relative to the pump
z1 = 2.0  # Height of water in the tank in m
z2 = 3.0  # Height of the pipe outlet in m

# Pipe properties
r_pipe = 15.5 / 1000  # Pipe radius in m (converted from mm)
L_pipe = 14.9  # Pipe length in m
f_friction = 0.004  # Darcy friction factor

# Minor loss coefficients
K_entrance = 0.4  # ef (shrinkage)

# --- 2. Calculate Intermediate Values ---
# Pipe diameter and area
D_pipe = 2 * r_pipe
A_pipe = math.pi * r_pipe**2

# Fluid velocity in the pipe
v_pipe = Q / A_pipe

# Velocity head
velocity_head = v_pipe**2 / (2 * g)

# --- 3. Calculate Head Losses ---
# Major loss (pipe friction)
# The L/D from the table (31) is inconsistent with L and r. We calculate it.
L_over_D = L_pipe / D_pipe
h_loss_major = f_friction * L_over_D * velocity_head

# Minor loss (entrance)
h_loss_minor_entrance = K_entrance * velocity_head

# Head loss in the cooling mold, calculated from the pressure drop
# h_loss = (P_in - P_out)/rho*g + (v_in^2 - v_out^2)/2g
# v_in = v_pipe, v_out is assumed negligible
h_loss_mold = (P_manometer - P_ambient) / (rho * g) + velocity_head

# Total head loss is the sum of all losses
h_loss_total = h_loss_major + h_loss_minor_entrance + h_loss_mold

# --- 4. Calculate Pump Head ---
# From Bernoulli equation: h_pump = (z2-z1) + h_total_loss
h_pump = (z2 - z1) + h_loss_total

# --- 5. Calculate Pump Work (Power) ---
# Power = rho * Q * g * h_pump
pump_work = rho * Q * g * h_pump

# --- 6. Print the Results ---
print("--- Calculation Steps ---")
print(f"Pipe velocity (v): {v_pipe:.3f} m/s")
print(f"Velocity head (v^2/2g): {velocity_head:.3f} m")
print(f"Static lift (z2 - z1): {z2 - z1:.3f} m")
print("\n--- Head Losses ---")
print(f"Major friction loss (pipe): {h_loss_major:.3f} m")
print(f"Minor loss (entrance): {h_loss_minor_entrance:.3f} m")
print(f"Head loss (cooling mold): {h_loss_mold:.3f} m")
print(f"Total head loss: {h_loss_total:.3f} m")
print("\n--- Pump Requirements ---")
print(f"Total pump head required (h_pump): {h_pump:.3f} m")

print("\n--- Final Equation for Pump Work (Power) ---")
print("Work_pump = ρ * Q * g * h_pump")
print("Work_pump = ρ * Q * g * [(z₂ - z₁) + h_L_major + h_L_minor_entrance + h_L_mold]")
print(f"Work_pump = {rho} * {Q} * {g} * [({z2} - {z1}) + {h_loss_major:.3f} + {h_loss_minor_entrance:.3f} + {h_loss_mold:.3f}]")
print(f"Work_pump = {rho} * {Q} * {g} * [{h_pump:.3f}]")
print(f"Work_pump = {pump_work:.2f} W")
print(f"\nThe calculated work of the pump is {pump_work:.2f} Watts.")

<<<175.91>>>