import math

# --- 1. Define known values from the problem statement ---
rho = 997  # Density of water in kg/m^3
g = 9.81  # Acceleration due to gravity in m/s^2
P1 = 1.0e5  # Ambient pressure at tank surface in N/m^2
PM = 1.33e5  # Manometer absolute pressure in N/m^2
Q = 2.86e-3  # Volume flow rate in m^3/s

# Geometric and system parameters
r_mm = 15.5  # Pipe radius in mm
h_water_tank = 2  # Height of water in the tank in m
h_pipe_vertical = 3  # Vertical length of the pipe in m
L_pipe = 14.9  # Total length of the pipe in m
Le_over_D_bend = 31  # Equivalent length ratio for the bend
f = 0.004  # Darcy friction factor
K_entrance = 0.4  # Minor loss coefficient for shrinkage (tank exit)

# --- 2. Convert units and calculate derived geometric properties ---
r = r_mm / 1000  # Convert radius to meters
D = 2 * r  # Pipe diameter in m
A = math.pi * r**2  # Pipe cross-sectional area in m^2

# --- 3. Calculate fluid velocity and related head terms ---
v_pipe = Q / A  # Velocity of water in the pipe in m/s
velocity_head = v_pipe**2 / (2 * g)  # Velocity head (v^2 / 2g) in m

# Define elevations
# Let's set the manometer elevation (zM) as the reference datum (0 m).
zM = 0
# The water surface (z1) is the sum of the vertical pipe length and water height.
z1 = h_water_tank + h_pipe_vertical
elevation_head_diff = z1 - zM # Elevation difference in m

# Calculate pressure head difference
pressure_head_diff = (PM - P1) / (rho * g)

# --- 4. Calculate total head loss from tank to manometer (h_friction) ---
# Major loss coefficient
major_loss_coeff = f * (L_pipe / D)
# Minor loss for the bend
K_bend = f * Le_over_D_bend
# Total head loss
h_friction = (major_loss_coeff + K_entrance + K_bend) * velocity_head

# --- 5. Rearrange Bernoulli's equation to solve for pump head (h_pump) ---
# h_pump = (P_M - P₁)/(ρg) + (v_M²/2g) - (z₁ - z_M) + h_friction
# We are adding energy with the pump, and overcoming pressure difference, elevation change, and friction.
h_pump = pressure_head_diff + velocity_head - elevation_head_diff + h_friction

# --- 6. Calculate the work of the pump (Power in Watts) ---
W_pump = rho * g * Q * h_pump

# --- 7. Print the results step-by-step ---
print("--- Calculation Steps ---")
print(f"1. Pump Head (h_pump) Calculation:")
print(f"h_pump = (Pressure Head Change) + (Velocity Head) - (Elevation Head Change) + (Friction Head Loss)")
print(f"h_pump = (({PM:.2e} - {P1:.2e}) / ({rho} * {g})) + ({v_pipe**2:.2f} / (2 * {g})) - ({z1} - {zM}) + {h_friction:.4f}")
print(f"h_pump = {pressure_head_diff:.4f} m + {velocity_head:.4f} m - {elevation_head_diff:.4f} m + {h_friction:.4f} m")
print(f"h_pump = {h_pump:.4f} m\n")

print(f"2. Pump Work (Power) Calculation:")
print(f"Work = ρ * g * Q * h_pump")
print(f"Work = {rho} kg/m^3 * {g} m/s^2 * {Q} m^3/s * {h_pump:.4f} m")
print(f"Work = {W_pump:.2f} Watts")

print("\n--- Final Answer ---")
print(f"The calculated work of the pump is {W_pump:.2f} Watts.")
<<<25.08>>>