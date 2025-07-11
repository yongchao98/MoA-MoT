import math

# 1. Define all the given variables from the problem in SI units.
rho = 997  # Density of water in kg/m^3
g = 9.81   # Acceleration due to gravity in m/s^2
Q = 2.86e-3 # Volume flow rate in m^3/s

# Pressures P1 and P2 are atmospheric and equal, so they cancel out.
P_ambient = 1.0e5 # N/m^2 (Assuming the typo 105 meant 1.0e5, value not needed as P1=P2)
P_manometer = 1.33e5 # N/m^2 (This data is considered secondary to friction parameters)

# Elevations relative to the pump level (z=0)
z1 = 2.0    # Height of water surface in the tank in m
z2 = -3.0   # Height of the exit in m

# Pipe and friction parameters from the table
r = 15.5 / 1000 # Pipe radius in m
L = 14.9        # Pipe length in m
L_D_bend = 31   # Equivalent L/D for the bend
f = 0.004       # Fanning friction factor
K_contraction = 0.4 # Loss coefficient for entrance (shrinkage)
K_expansion = 0.8   # Loss coefficient for exit (expansion)

# 2. Calculate intermediate values
D = 2 * r  # Pipe diameter in m
A = math.pi * r**2 # Pipe cross-sectional area in m^2
v2 = Q / A  # Fluid velocity at exit (and in the pipe) in m/s
v2_squared_over_2g = v2**2 / (2 * g) # Kinetic energy head in m

# 3. Calculate total head loss (h_L)
# Head loss from straight pipe friction
hL_pipe_factor = 4 * f * (L / D)
hL_pipe = hL_pipe_factor * v2_squared_over_2g

# Head loss from the bend, using equivalent L/D
hL_bend_factor = 4 * f * L_D_bend
hL_bend = hL_bend_factor * v2_squared_over_2g

# Head loss from fittings (entrance contraction and exit expansion)
hL_fittings_factor = K_contraction + K_expansion
hL_fittings = hL_fittings_factor * v2_squared_over_2g

# Total head loss
hL_total = hL_pipe + hL_bend + hL_fittings

# 4. Solve for the pump head (h_p) using the rearranged Bernoulli equation
# h_p = (z2 - z1) + (v2^2 / 2g) + h_L
delta_z = z2 - z1
h_p = delta_z + v2_squared_over_2g + hL_total

# 5. Calculate the work of the pump (Power) in Watts
Work_pump = rho * g * Q * h_p

# --- Output the results and the final equations ---
print("--- Step-by-step Calculation ---")
print(f"Pipe Diameter (D): {D:.4f} m")
print(f"Pipe Area (A): {A:.6f} m^2")
print(f"Fluid Velocity (v2): {v2:.3f} m/s")
print(f"Kinetic Head (v2^2/2g): {v2_squared_over_2g:.3f} m")
print("\n--- Head Loss Calculation ---")
print(f"Pipe Friction Loss (hL_pipe): {hL_pipe:.3f} m")
print(f"Bend Friction Loss (hL_bend): {hL_bend:.3f} m")
print(f"Fittings Loss (hL_fittings): {hL_fittings:.3f} m")
print(f"Total Head Loss (hL_total): {hL_total:.3f} m")
print("\n--- Pump Head Calculation ---")
print("The pump head (hp) is found from the Bernoulli equation: hp = (z2 - z1) + v2^2/(2g) + hL_total")
print(f"hp = ({z2:.1f} m - {z1:.1f} m) + {v2_squared_over_2g:.3f} m + {hL_total:.3f} m")
print(f"hp = {delta_z:.1f} m + {v2_squared_over_2g:.3f} m + {hL_total:.3f} m = {h_p:.3f} m")
print("\n--- Pump Work (Power) Calculation ---")
print("The work of the pump (Power) is calculated as: Power = ρ * g * Q * hp")
print(f"Power = {rho} kg/m^3 * {g} m/s^2 * {Q} m^3/s * {h_p:.3f} m")
print(f"Power = {Work_pump:.2f} Watts")

# Final answer in the required format
final_answer = Work_pump
print(f"\nFinal calculated work of the pump is {final_answer:.2f} Watts.")
print(f"\n<<<72.72>>>")