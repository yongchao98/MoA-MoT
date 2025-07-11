import math

# --- 1. Define constants and given values ---
rho = 997      # Density of water (kg/m^3)
P_amb = 1e5    # Ambient pressure (N/m^2), assuming 105 is a typo for 10^5
P_p = 1.33e5   # Absolute pressure at manometer P (N/m^2)
Q = 2.86e-3    # Volume flow rate (m^3/s)
g = 9.81       # Acceleration due to gravity (m/s^2)
delta_z = -3   # Elevation change (z2 - z1) in meters

# Pipe and friction data from the table
r = 15.5e-3    # Pipe radius (m)
L_pipe = 14.9  # Pipe length (m)
L_D_fitting = 31 # L/D for fittings
f_fanning = 0.004 # Fanning friction factor
K_entrance = 0.4  # Loss coefficient for shrinkage (entrance)

# --- 2. Intermediate Calculations ---
D = 2 * r      # Pipe diameter (m)
A = math.pi * r**2 # Pipe cross-sectional area (m^2)
v = Q / A      # Fluid velocity in the pipe (m/s)
v_head = v**2 / (2 * g) # Velocity head (m)

# --- 3. Calculate Head Loss before Manometer (h_L_1_P) ---
L_D_pipe = L_pipe / D
K_major = 4 * f_fanning * (L_D_pipe + L_D_fitting)
K_total_1_P = K_major + K_entrance
h_L_1_P = K_total_1_P * v_head

# --- 4. Calculate Head Loss after Manometer (h_L_P_2) ---
# Apply energy equation from P to 2: (Pp/ρg) + (vp²/2g) = (P2/ρg) + (v2²/2g) + h_L(P->2)
# With v2=0, h_L(P->2) = (Pp - P2)/ρg + vp²/2g
pressure_head_diff = (P_p - P_amb) / (rho * g)
h_L_P_2 = pressure_head_diff + v_head

# --- 5. Calculate Total Head Loss, Pump Head, and Pump Work ---
h_L_total = h_L_1_P + h_L_P_2
h_pump = h_L_total + delta_z
work_pump = g * h_pump

# --- 6. Print the results and the final equation ---
print("--- Calculation Steps ---")
print(f"Pipe velocity (v): {v:.3f} m/s")
print(f"Velocity head (v^2/2g): {v_head:.3f} m")
print(f"Head loss before manometer (h_L_1_P): {h_L_1_P:.3f} m")
print(f"Head loss in mold (h_L_P_2): {h_L_P_2:.3f} m")
print(f"Total head loss (h_L): {h_L_total:.3f} m")
print(f"Elevation change (z2 - z1): {delta_z:.3f} m")
print(f"Required pump head (h_pump): {h_pump:.3f} m")
print("\n--- Final Work Calculation ---")
print("The work of the pump (W_pump) is calculated as: W_pump = g * (h_L_total + delta_z)")
print(f"W_pump = {g:.2f} m/s^2 * ({h_L_total:.3f} m + ({delta_z:.3f} m))")
print(f"W_pump = {g:.2f} * {h_pump:.3f}")
print(f"Work of the pump = {work_pump:.2f} J/kg")

print(f"\nFinal Answer:")
print(f"<<<{work_pump:.2f}>>>")