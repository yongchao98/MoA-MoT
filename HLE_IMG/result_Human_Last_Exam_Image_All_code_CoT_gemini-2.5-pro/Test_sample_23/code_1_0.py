import math

# --- Given Data ---
rho = 997  # Density of water in kg/m^3
g = 9.81   # Acceleration due to gravity in m/s^2
P_ambient = 1.05e5  # Ambient pressure in N/m^2 (assuming 105 kPa from the given 105)
P_manometer = 1.33e5 # Absolute pressure at manometer in N/m^2
Q = 2.86e-3  # Volume flow rate in m^3/s
h_tank = 2   # Height of water in the tank above the pump in m
r_pipe = 15.5 / 1000 # Pipe radius in m (15.5 mm)
L_over_D = 31 # Length to diameter ratio of the pipe
f = 0.004  # Friction factor
ef_shrinkage = 0.4 # Minor loss coefficient for pipe entrance (shrinkage)

# --- Plan ---
# 1. Define points A (tank surface) and M (manometer).
# 2. Use Bernoulli's equation: (Pa/rho*g) + (va^2/2g) + za + H_pump = (Pm/rho*g) + (vm^2/2g) + zm + H_loss
# 3. Rearrange to solve for H_pump: H_pump = (Pm-Pa)/(rho*g) + (vm^2/2g) + (zm-za) + H_loss
# 4. Calculate required values: vm (pipe velocity) and H_loss.
# 5. Calculate H_pump (pump head).
# 6. Calculate the work of the pump (Power) = rho * g * Q * H_pump.

# --- Calculations ---

# Pipe cross-sectional area (A)
A_pipe = math.pi * r_pipe**2

# Velocity in the pipe (vm)
v_pipe = Q / A_pipe

# Elevation head difference (zm - za)
# Let pump level (zm) be the datum (0 m). The tank surface (za) is 2 m above.
z_diff = 0 - h_tank

# Pressure head difference
pressure_head = (P_manometer - P_ambient) / (rho * g)

# Velocity head at the manometer
velocity_head = (v_pipe**2) / (2 * g)

# Head loss from tank to manometer (H_loss)
# H_loss = (Loss due to entrance) + (Loss due to friction)
# H_loss = K_entrance * (v^2/2g) + f*(L/D)*(v^2/2g) = (K_entrance + f*(L/D)) * (v^2/2g)
head_loss = (ef_shrinkage + f * L_over_D) * velocity_head

# Pump Head (H_pump)
H_pump = pressure_head + velocity_head + z_diff + head_loss

# Work of the pump (Power in Watts)
W_pump_power = rho * g * Q * H_pump

# --- Output Results ---
print("--- Calculation of Pump Head (H_pump) ---")
print("The pump head is calculated using the Bernoulli equation rearranged as:")
print("H_pump = (Pressure Head) + (Velocity Head) + (Elevation Head) + (Head Loss)\n")

print(f"Pressure Head = (P_manometer - P_ambient) / (rho * g)")
print(f"Pressure Head = ({P_manometer:.0f} - {P_ambient:.0f}) N/m^2 / ({rho} kg/m^3 * {g} m/s^2) = {pressure_head:.2f} m")

print(f"\nVelocity Head = v_pipe^2 / (2 * g)")
print(f"Velocity Head = {v_pipe:.2f}^2 m^2/s^2 / (2 * {g} m/s^2) = {velocity_head:.2f} m")

print(f"\nElevation Head = z_manometer - z_tank_surface = {z_diff:.2f} m")

print(f"\nHead Loss = (ef_shrinkage + f * L/D) * Velocity Head")
print(f"Head Loss = ({ef_shrinkage} + {f} * {L_over_D}) * {velocity_head:.2f} m = {head_loss:.2f} m")

print("\n--- Total Pump Head ---")
print(f"H_pump = {pressure_head:.2f} m + {velocity_head:.2f} m + ({z_diff:.2f} m) + {head_loss:.2f} m")
print(f"H_pump = {H_pump:.2f} m\n")


print("--- Final Calculation of the Work of the Pump (Power) ---")
print("Power = rho * g * Q * H_pump")
print(f"Power = {rho} kg/m^3 * {g} m/s^2 * {Q} m^3/s * {H_pump:.2f} m")
print(f"Power = {W_pump_power:.2f} Watts")

# Final answer in the requested format
# The final answer is the work of the pump (Power) in Watts.
final_answer = W_pump_power
#<<<55.33>>>