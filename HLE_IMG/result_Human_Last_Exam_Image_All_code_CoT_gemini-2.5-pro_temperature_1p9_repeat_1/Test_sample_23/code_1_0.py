import math

# --- 1. Define constants and given parameters ---
rho = 997  # Density of water in kg/m^3
g = 9.81  # Acceleration due to gravity in m/s^2
P1 = 1.0e5  # Pressure at tank surface (ambient) in N/m^2
Pm = 1.33e5 # Pressure at manometer in N/m^2
Q = 2.86e-3 # Volume flow rate in m^3/s
r_mm = 15.5 # Pipe radius in mm
L_pipe = 14.9 # Pipe length in m
f_fanning = 0.004 # Fanning friction factor
LD_eq_fittings = 31 # L/D equivalent for fittings
K_shrinkage = 0.4 # Loss coefficient for entrance (shrinkage)
z1 = 5 # Elevation of water surface in tank in m (2m water + 3m pipe drop)
zm = 0 # Elevation of manometer (datum) in m
v1 = 0 # Velocity at water surface (negligible)

# --- 2. Calculate intermediate values ---
# Convert radius to meters and calculate diameter
r_m = r_mm / 1000.0
D = 2 * r_m

# Calculate Darcy friction factor
f_darcy = 4 * f_fanning

# Calculate pipe cross-sectional area
A = math.pi * r_m**2

# Calculate fluid velocity in the pipe (vm)
vm = Q / A

# --- 3. Calculate head loss (h_L) from tank to manometer ---
# Head loss is due to entrance, pipe friction, and fittings
# Loss coefficient for fittings
K_fittings = f_darcy * LD_eq_fittings

# Total head loss term calculation
velocity_head_term = (vm**2) / (2 * g)
major_loss = f_darcy * (L_pipe / D) * velocity_head_term
minor_loss_entrance = K_shrinkage * velocity_head_term
minor_loss_fittings = K_fittings * velocity_head_term

h_L = major_loss + minor_loss_entrance + minor_loss_fittings

# --- 4. Apply Bernoulli's equation to find pump head (h_p) ---
pressure_head = (Pm - P1) / (rho * g)
velocity_head = (vm**2 - v1**2) / (2 * g)
elevation_head = zm - z1

h_p = pressure_head + velocity_head + elevation_head + h_L

# --- 5. Calculate pump power (Work Rate) in Watts ---
# The "work of the pump" is typically interpreted as power
pump_power = rho * g * Q * h_p

# --- 6. Print the results step-by-step ---
print("--- Calculation Steps ---")
print(f"1. Fluid velocity in pipe (v_m): {vm:.2f} m/s")
print(f"2. Darcy friction factor (f_D): {f_darcy:.4f}")
print(f"3. Head Loss Calculation:")
print(f"   - Major loss (pipe friction): {major_loss:.2f} m")
print(f"   - Minor loss (entrance): {minor_loss_entrance:.2f} m")
print(f"   - Minor loss (fittings): {minor_loss_fittings:.2f} m")
print(f"   - Total Head Loss (h_L): {h_L:.2f} m")
print(f"4. Pump Head Calculation (h_p):")
print(f"   - Pressure head term: {pressure_head:.2f} m")
print(f"   - Velocity head term: {velocity_head:.2f} m")
print(f"   - Elevation head term: {elevation_head:.2f} m")
print(f"   h_p = {pressure_head:.2f} + {velocity_head:.2f} + ({elevation_head:.2f}) + {h_L:.2f} = {h_p:.2f} m")
print("\n--- Final Answer ---")
print("The final equation for pump power (Work Rate) is:")
print(f"Work_pump = rho * g * Q * h_p")
print(f"Work_pump = {rho} kg/m^3 * {g} m/s^2 * {Q} m^3/s * {h_p:.2f} m")
print(f"Work of the pump = {pump_power:.2f} Watts")

<<<150.78>>>