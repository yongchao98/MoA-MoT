import math

# --- Given Parameters ---
rho = 997  # Density of water in kg/m^3
g = 9.81   # Acceleration due to gravity in m/s^2
Q = 2.86e-3 # Volume flow rate in m^3/s
z1 = 2     # Elevation of water in the tank in meters
z2 = 3     # Elevation of the pipe outlet in meters
r_pipe_mm = 15.5 # Pipe radius in mm
L_pipe = 14.9    # Pipe length in meters
f = 0.004        # Darcy friction factor
K_in = 0.4       # Minor loss coefficient for entrance (shrinkage)
K_out = 0.8      # Minor loss coefficient for exit (expansion)
P_atm = 1.0e5    # Ambient pressure in N/m^2
# The manometer pressure is noted but will not be used in the primary calculation
# as it conflicts with the theoretical loss calculation parameters provided.

# --- Step 1: Calculate pipe properties and fluid velocity ---
r_pipe = r_pipe_mm / 1000  # Convert radius to meters
D_pipe = 2 * r_pipe
A_pipe = math.pi * r_pipe**2
v_pipe = Q / A_pipe

print(f"Pipe diameter (D): {D_pipe:.4f} m")
print(f"Pipe area (A): {A_pipe:.6f} m^2")
print(f"Water velocity in pipe (v): {v_pipe:.4f} m/s")

# --- Step 2: Calculate the kinetic head term ---
v_sq_over_2g = v_pipe**2 / (2 * g)
print(f"Kinetic head (v^2/2g): {v_sq_over_2g:.4f} m")

# --- Step 3: Calculate the total head loss (H_L) ---
# Major loss (friction)
H_L_major = f * (L_pipe / D_pipe) * v_sq_over_2g

# Minor losses (fittings)
sum_K = K_in + K_out
H_L_minor = sum_K * v_sq_over_2g

# Total head loss
H_L_total = H_L_major + H_L_minor

print(f"\nHead loss from friction: {H_L_major:.4f} m")
print(f"Head loss from fittings: {H_L_minor:.4f} m")
print(f"Total head loss (H_L): {H_L_total:.4f} m")

# --- Step 4: Calculate the required pump head (H_p) ---
delta_z = z2 - z1
H_pump = delta_z + H_L_total

print(f"\nRequired elevation gain (z2 - z1): {delta_z:.4f} m")
print(f"Required pump head (H_p = (z2-z1) + H_L): {H_pump:.4f} m")

# --- Step 5: Calculate the work (power) of the pump ---
Work_pump = rho * g * Q * H_pump

print("\n--- Final Calculation for Work of the Pump (Power) ---")
print(f"Work = ρ * g * Q * H_p")
print(f"Work = {rho} kg/m^3 * {g} m/s^2 * {Q} m^3/s * {H_pump:.4f} m")
print(f"Work = {rho * g * Q:.4f} N/s * {H_pump:.4f} m")
print(f"Work = {Work_pump:.4f} Watts")

<<<91.9023>>>