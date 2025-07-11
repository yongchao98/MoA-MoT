import math

# --- 1. Define Given Variables and Constants ---
rho = 997  # Density of water in kg/m^3
g = 9.81  # Acceleration of gravity in m/s^2
p_ambient = 1.0e5  # Ambient pressure at tank surface and exit in N/m^2
p_manometer = 1.33e5  # Absolute pressure at the manometer in N/m^2
Q = 2.86e-3  # Volume flow rate in m^3/s
r_pipe = 15.5 / 1000  # Pipe radius in m
L_over_D = 31  # L/D ratio for the pipe
f = 0.004  # Darcy friction factor for the pipe
k_entrance = 0.4  # Entrance loss coefficient (e_f shrinkage)
delta_z_1_to_2 = -3.0 # Elevation change from tank surface to exit (z2 - z1) in m

print("This script calculates the work of the pump (Power) based on the provided system parameters.")
print("The calculation is performed step-by-step.\n")

# --- 2. Calculate Fluid Velocity and Velocity Head ---
d_pipe = 2 * r_pipe
A_pipe = math.pi * r_pipe**2
v_pipe = Q / A_pipe
v_head = (v_pipe**2) / (2 * g)

print("--- Initial Calculations ---")
print(f"Pipe Diameter (D): {d_pipe:.4f} m")
print(f"Pipe Area (A): {A_pipe:.6f} m^2")
print(f"Water Velocity in Pipe (v): {v_pipe:.4f} m/s")
print(f"Velocity Head (v^2 / 2g): {v_head:.4f} m\n")

# --- 3. Calculate Head Loss in the Cooling Mold ---
# Applying Bernoulli from manometer (m) to exit (2), assuming zm=z2 and v2=0
# H_loss_mold = (Pm - P2)/rho*g + (vm^2 - v2^2)/2g + (zm-z2)
h_loss_mold = (p_manometer - p_ambient) / (rho * g) + v_head
print("--- Step 1: Calculate Head Loss in the Mold ---")
print("Equation: H_loss_mold = (P_manometer - P_ambient) / (rho * g) + (v_pipe^2 / (2 * g))")
print(f"H_loss_mold = ({p_manometer:.0f} - {p_ambient:.0f}) / ({rho} * {g}) + {v_head:.4f}")
print(f"H_loss_mold = {h_loss_mold:.4f} m\n")

# --- 4. Calculate Other Head Losses (Entrance and Pipe) ---
h_loss_entrance = k_entrance * v_head
h_loss_pipe = f * L_over_D * v_head
print("--- Step 2: Calculate Other System Head Losses ---")
print("Equation: H_loss_entrance = k_entrance * (v_pipe^2 / (2 * g))")
print(f"H_loss_entrance = {k_entrance} * {v_head:.4f} = {h_loss_entrance:.4f} m")
print("\nEquation: H_loss_pipe = f * (L/D) * (v_pipe^2 / (2 * g))")
print(f"H_loss_pipe = {f} * {L_over_D} * {v_head:.4f} = {h_loss_pipe:.4f} m\n")


# --- 5. Calculate Total Head Loss ---
h_loss_total = h_loss_entrance + h_loss_pipe + h_loss_mold
print("--- Step 3: Calculate Total Head Loss ---")
print("Equation: H_loss_total = H_loss_entrance + H_loss_pipe + H_loss_mold")
print(f"H_loss_total = {h_loss_entrance:.4f} m + {h_loss_pipe:.4f} m + {h_loss_mold:.4f} m")
print(f"H_loss_total = {h_loss_total:.4f} m\n")

# --- 6. Calculate Pump Head ---
# Applying overall Bernoulli from tank surface (1) to exit (2)
# P1/rhog + v1^2/2g + z1 + H_pump = P2/rhog + v2^2/2g + z2 + H_loss_total
# With P1=P2, v1=0, v2=0, this simplifies to:
# H_pump = (z2 - z1) + H_loss_total
H_pump = delta_z_1_to_2 + h_loss_total
print("--- Step 4: Calculate Pump Head (H_pump) ---")
print("Equation: H_pump = (z_exit - z_tank) + H_loss_total")
print(f"H_pump = {delta_z_1_to_2:.4f} m + {h_loss_total:.4f} m")
print(f"H_pump = {H_pump:.4f} m\n")

# --- 7. Calculate Pump Work (Power) ---
W_pump = H_pump * rho * g * Q
print("--- Step 5: Calculate the Work of the Pump (Power) ---")
print("Equation: Power = H_pump * rho * g * Q")
print(f"Power = {H_pump:.4f} m * {rho} kg/m^3 * {g} m/s^2 * {Q} m^3/s")
print(f"\nThe work of the pump is {W_pump:.2f} Watts.")

# The final answer as a number
final_answer = W_pump
# print(f"\n<<<{final_answer:.2f}>>>")