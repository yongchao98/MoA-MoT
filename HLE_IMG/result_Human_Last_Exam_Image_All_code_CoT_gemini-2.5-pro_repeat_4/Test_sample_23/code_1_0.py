import math

# Step 1: Define given data and constants
rho = 997  # Density of water in kg/m^3
P_ambient = 1.05e5  # Ambient pressure in N/m^2
Q = 2.86e-3  # Volume flow rate in m^3/s
z1_relative = 2  # Height of water in tank above the pump in m
z2_relative = 3  # Vertical length of the pipe from the pump to the exit in m
r_mm = 15.5  # Pipe radius in mm
L_D_ratio = 31  # Given L/D ratio for friction calculation
f = 0.004  # Darcy friction factor
K_shrinkage = 0.4  # Minor loss coefficient for entrance (shrinkage)
K_expansion = 0.8  # Minor loss coefficient for exit/mold (expansion)
g = 9.81  # Acceleration due to gravity in m/s^2

# Convert radius to meters
r_m = r_mm / 1000

print("--- Step-by-step Calculation for the Work of the Pump ---")

# Step 2: Calculate pipe cross-sectional area (A)
A = math.pi * r_m**2
print(f"\n1. Pipe Cross-Sectional Area (A):")
print(f"A = π * r^2 = π * ({r_m})^2 = {A:.6f} m^2")

# Step 3: Calculate water velocity in the pipe (v)
v = Q / A
print(f"\n2. Water Velocity (v):")
print(f"v = Q / A = {Q} / {A:.6f} = {v:.4f} m/s")

# Step 4: Calculate the velocity head (v^2 / 2g)
velocity_head = v**2 / (2 * g)
print(f"\n3. Velocity Head (v^2 / 2g):")
print(f"v^2 / (2*g) = ({v:.4f})^2 / (2 * {g}) = {velocity_head:.4f} m")

# Step 5: Calculate the total head loss (h_L)
# h_L = Major Loss (friction) + Minor Losses (fittings)
h_friction = f * L_D_ratio * velocity_head
h_minor = (K_shrinkage + K_expansion) * velocity_head
h_L = h_friction + h_minor
print(f"\n4. Total Head Loss (h_L):")
print(f"h_L = (f * (L/D) + K_shrinkage + K_expansion) * (v^2 / 2g)")
print(f"h_L = ({f} * {L_D_ratio} + {K_shrinkage} + {K_expansion}) * {velocity_head:.4f}")
print(f"h_L = ({f*L_D_ratio:.3f} + {K_shrinkage + K_expansion}) * {velocity_head:.4f}")
print(f"h_L = {f*L_D_ratio + K_shrinkage + K_expansion:.3f} * {velocity_head:.4f} = {h_L:.4f} m")

# Step 6: Calculate the pump head (H_p) using the simplified Bernoulli equation
# H_p = (P2-P1)/ρg + (v2^2-v1^2)/2g + (z2-z1) + h_L
# With P1=P2 and v1=v2=0, this simplifies to H_p = (z2-z1) + h_L
delta_z = z2_relative - z1_relative
H_p = delta_z + h_L
print(f"\n5. Pump Head (H_p):")
print(f"H_p = (z2 - z1) + h_L")
print(f"H_p = ({z2_relative} - {z1_relative}) + {h_L:.4f}")
print(f"H_p = {delta_z} + {h_L:.4f} = {H_p:.4f} m")

# Step 7: Calculate the specific work of the pump (W_pump)
W_pump = H_p * g
print(f"\n6. Specific Work of the Pump (W_pump):")
print(f"W_pump = H_p * g")
print(f"W_pump = {H_p:.4f} * {g} = {W_pump:.2f} J/kg")
print(f"\nThe calculated work of the pump is {W_pump:.2f} J/kg.")
