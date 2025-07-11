import math

# --- Given Data ---
rho = 997  # Density of water in kg/m^3
P1 = 1.05e5  # Ambient pressure at source tank surface in N/m^2
P_m = 1.33e5 # Manometer absolute pressure in N/m^2
Q = 2.86e-3  # Volume flow rate in m^3/s
z1 = 2.0  # Height of water in the tank relative to the pump in m
z_m = 3.0 # Vertical length of the pipe (outlet height) relative to the pump in m
g = 9.81 # Acceleration due to gravity in m/s^2
r = 15.5 / 1000 # Pipe radius in m
L = 14.9 # Pipe length in m
f = 0.004 # Darcy friction factor
K_c = 0.4 # Minor loss coefficient for entrance (shrinkage)
K_e = 0.8 # Minor loss coefficient for exit (expansion)

# Note on data consistency:
# There is a conflict between the pressure drop calculated from the manometer reading 
# (P_m - P_ambient) and the pressure drop calculated using the exit loss coefficient K_e.
# We will proceed using the manometer reading (P_m) as it represents a direct measurement
# for this specific system. The Bernoulli equation will be applied from the tank surface (1)
# to the manometer location (m).

# --- Calculations ---

# 1. Calculate pipe diameter and cross-sectional area
D = 2 * r
A = math.pi * r**2

# 2. Calculate water velocity in the pipe
v_m = Q / A

# 3. Calculate terms for the Bernoulli equation

# Pressure head term
pressure_term = (P_m - P1) / rho

# Kinetic energy term (since v1 is negligible)
kinetic_term = (v_m**2) / 2

# Elevation head term
elevation_term = g * (z_m - z1)

# 4. Calculate total head loss from tank to manometer
# Major loss (pipe friction)
major_loss = f * (L / D) * ((v_m**2) / 2)
# Minor loss (entrance)
minor_loss = K_c * ((v_m**2) / 2)
# Total loss to manometer
h_L_total = major_loss + minor_loss

# 5. Calculate work done by the pump per unit mass (w_p)
# w_p = (P_m - P1)/rho + v_m²/2 + g*(z_m - z1) + h_L_total
w_p = pressure_term + kinetic_term + elevation_term + h_L_total

# 6. Calculate the total power of the pump (Work per unit time)
pump_power = w_p * rho * Q

# --- Output the Results ---
print("--- Calculation Steps ---")
print(f"1. Pipe velocity (v_m):")
print(f"   v_m = Q / (pi*r^2) = {Q:.3g} / (3.14159 * {r:.3g}^2) = {v_m:.3f} m/s\n")

print("2. Work per unit mass (w_p) using Bernoulli's equation:")
print("   w_p = (P_m - P1)/rho + v_m^2/2 + g*(z_m - z1) + [f*(L/D) + K_c]*v_m^2/2")
print(f"   w_p = ({P_m:.2e} - {P1:.2e})/{rho} + {v_m:.2f}^2/2 + {g:.2f}*({z_m} - {z1}) + [{f}*({L}/{D:.3f}) + {K_c}]*{v_m:.2f}^2/2")
print(f"   w_p = {pressure_term:.2f} + {kinetic_term:.2f} + {elevation_term:.2f} + [{major_loss:.2f} + {minor_loss:.2f}]")
print(f"   w_p = {pressure_term:.2f} + {kinetic_term:.2f} + {elevation_term:.2f} + {h_L_total:.2f}")
print(f"   w_p = {w_p:.2f} J/kg\n")

print("3. Total work of the pump (Power):")
print("   Power = w_p * rho * Q")
print(f"   Power = {w_p:.2f} * {rho} * {Q:.3g}")
print(f"   Power = {pump_power:.2f} W")

# Final answer in the required format
final_answer = round(pump_power)
print(f"\nThe calculated work of the pump is {final_answer} W.")
<<<176>>>