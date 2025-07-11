import math

# --- 1. Define given variables and constants ---
rho = 997  # Density of water in kg/m^3
g = 9.81  # Acceleration due to gravity in m/s^2
Q = 2.86e-3  # Volume flow rate in m^3/s

# Pipe properties
r_pipe = 15.5 / 1000  # Pipe radius in meters
L_pipe = 14.9  # Pipe length in meters
D_pipe = 2 * r_pipe

# Coefficients from the table
# Assumption: The provided friction factor f=0.004 is a Fanning friction factor.
# The Darcy-Weisbach equation uses the Darcy friction factor (f_D = 4 * f_Fanning).
f_fanning = 0.004
f_darcy = 4 * f_fanning
K_contraction = 0.4  # ef (shrinkage)
K_expansion = 0.8  # ef (expansion)

# System elevations
# Assumption: The initial height (z1) is the sum of water in the tank and the vertical pipe length.
# The final height (z2) is the reference datum (0).
h_water_tank = 2.0  # meters
h_vertical_pipe = 3.0  # meters
z1 = h_water_tank + h_vertical_pipe
z2 = 0
delta_z = z2 - z1

# --- 2. Calculate intermediate values ---
# Pipe cross-sectional area
A_pipe = math.pi * r_pipe**2

# Water velocity in the pipe
v_pipe = Q / A_pipe

# Velocity head
velocity_head = v_pipe**2 / (2 * g)

# --- 3. Calculate head losses ---
# Major head loss (due to friction in the straight pipe)
# We calculate L/D from the given L and r, ignoring the inconsistent table value.
L_over_D = L_pipe / D_pipe
h_f_major = f_darcy * L_over_D * velocity_head

# Minor head loss (due to entrance contraction)
h_f_contraction = K_contraction * velocity_head

# Minor head loss (due to exit expansion)
h_f_expansion = K_expansion * velocity_head

# Total head loss
h_f_total = h_f_major + h_f_contraction + h_f_expansion

# --- 4. Apply Bernoulli's equation to find pump head (h_p) ---
# The extended Bernoulli equation is:
# (P1/ρg) + (v1²/2g) + z1 + h_p = (P2/ρg) + (v2²/2g) + z2 + h_f_total
# Assuming P1=P2 and v1=v2=0, it simplifies to: z1 + h_p = z2 + h_f_total
# Rearranging for pump head: h_p = (z2 - z1) + h_f_total
h_p = delta_z + h_f_total

# --- 5. Calculate the work of the pump per unit mass (w_p) ---
w_p = h_p * g

# --- 6. Print the results step-by-step ---
print("--- Calculation Steps ---")
print(f"Pipe diameter (D): {D_pipe:.4f} m")
print(f"Pipe velocity (v): {v_pipe:.4f} m/s")
print(f"Velocity head (v^2 / 2g): {velocity_head:.4f} m")
print("\n--- Head Loss Calculation ---")
print(f"Major loss (friction): {h_f_major:.4f} m")
print(f"Minor loss (contraction): {h_f_contraction:.4f} m")
print(f"Minor loss (expansion): {h_f_expansion:.4f} m")
print(f"Total head loss (h_f_total): {h_f_total:.4f} m")
print("\n--- Pump Head and Work Calculation ---")
print("The final equation for the work of the pump per unit mass (w_p) is derived from the pump head (h_p):")
print("h_p = (z2 - z1) + h_f_total")
print(f"h_p = ({z2} - {z1}) + {h_f_total:.4f} = {h_p:.4f} m")
print("\nw_p = h_p * g")
print(f"w_p = {h_p:.4f} m * {g} m/s^2")
print(f"\nThe work of the pump is: {w_p:.2f} J/kg")