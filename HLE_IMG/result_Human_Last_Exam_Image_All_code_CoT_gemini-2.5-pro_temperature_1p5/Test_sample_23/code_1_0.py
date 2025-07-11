import math

# --- 1. Define Variables and Constants from the problem statement ---
rho = 997  # Density of water in kg/m^3
g = 9.81  # Acceleration due to gravity in m/s^2
Q = 2.86e-3  # Volume flow rate in m^3/s

# Geometric and pipe properties
r_pipe = 15.5 / 1000  # Pipe radius in meters
L_pipe = 14.9  # Pipe length in meters
z1 = 2.0  # Initial height of water in meters
z2 = 3.0  # Final height of water at exit in meters

# Friction and loss coefficients
f_darcy = 0.004  # Darcy friction factor
Kc_entrance = 0.4  # Minor loss coefficient for entrance (shrinkage)
Ke_exit = 0.8  # Minor loss coefficient for exit (expansion)

# Pressures at start and end points are ambient, so P1=P2.

# --- 2. Calculate Intermediate Values ---
# Pipe diameter and area
D_pipe = 2 * r_pipe
A_pipe = math.pi * r_pipe**2

# Water velocity in the pipe (v2, since v1 in the large tank is ~0)
v2 = Q / A_pipe

# Kinetic head (velocity head) component of the energy equation
velocity_head = (v2**2) / (2 * g)

# --- 3. Calculate Head Losses (h_L) ---
# Major head loss due to pipe friction
L_over_D = L_pipe / D_pipe
h_L_major = f_darcy * L_over_D * velocity_head

# Minor head losses from fittings (entrance and exit)
h_L_minor = (Kc_entrance + Ke_exit) * velocity_head

# Total head loss
h_L_total = h_L_major + h_L_minor

# --- 4. Solve for Pump Head (H_pump) ---
# Elevation change
delta_z = z2 - z1

# Using the rearranged Bernoulli equation: H_pump = (z2-z1) + (v2^2/2g) + h_L
# The pressure term is zero as P1 = P2 (both ambient)
H_pump = delta_z + velocity_head + h_L_total

# --- 5. Calculate Pump Work (Power) ---
# Work of the pump (Power in Watts) is W = H_pump * rho * Q * g
W_pump = H_pump * rho * Q * g

# --- 6. Print the results step-by-step ---
print("--- Calculation Steps ---")
print(f"Pipe Velocity (v2): {v2:.2f} m/s")
print(f"Elevation Head (z2 - z1): {delta_z:.2f} m")
print(f"Velocity Head (v2^2 / 2g): {velocity_head:.2f} m")
print(f"Major Head Loss: {h_L_major:.2f} m")
print(f"Minor Head Loss: {h_L_minor:.2f} m")
print(f"Total Head Loss (h_L): {h_L_total:.2f} m")
print(f"Total Pump Head Required (H_pump): {H_pump:.2f} m")
print("\n--- Final Calculation for Pump Work (Power) ---")
print("Work_pump = (delta_z + velocity_head + h_L_total) * rho * Q * g")
print(f"Work_pump = ({delta_z:.2f} m + {velocity_head:.2f} m + {h_L_total:.2f} m) * {rho} kg/m^3 * {Q} m^3/s * {g} m/s^2")
print(f"Work_pump = ({H_pump:.2f} m) * {rho} kg/m^3 * {Q} m^3/s * {g} m/s^2")
print(f"\nThe calculated work of the pump is: {W_pump:.2f} Watts")
print(f"<<<{W_pump:.2f}>>>")