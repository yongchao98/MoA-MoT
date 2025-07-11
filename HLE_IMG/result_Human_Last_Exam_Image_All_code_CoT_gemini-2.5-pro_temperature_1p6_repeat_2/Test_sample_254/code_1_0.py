import math

# --- Given Data ---
# Power of the source in Watts
P = 1e9  # 1 GW
# Area of the satellite mirror in m^2 (not directly needed in the final formula, but part of the conceptual derivation)
A_mirror = 1e6 # 1 km^2
# Area of the photovoltaic cell in m^2
S = 10.0
# Orbital period in seconds
T = 12.0 * 3600.0
# Lunar mass in kg
M = 7.35e22
# Lunar radius in meters
R = 1738e3
# Gravitational constant in m^3 kg^-1 s^-2
G = 6.67e-11

# --- Step 1: Calculate the semi-major axis 'a' of the orbit ---
# Using Kepler's Third Law: T^2 = (4 * pi^2 * a^3) / (G * M)
# Rearranging for a: a = (G * M * T^2 / (4 * pi^2))^(1/3)
pi = math.pi
a_cubed = (G * M * T**2) / (4 * pi**2)
a = a_cubed**(1/3)

print(f"Calculated semi-major axis (a): {a:.4e} m")

# --- Step 2: Determine the distance r_YB ---
# Based on the physics of the reflection, the satellite Y is at a distance 'a' from the Moon's center.
# Point B is on the surface directly below Y.
# r_YB = a - R
r_YB = a - R
print(f"Distance from satellite Y to point B (r_YB): {r_YB:.4e} m")

# --- Step 3: Calculate the incident power P' ---
# The incidence angle alpha_B is 0, so cos(alpha_B) = 1.
# The formula for the power P' is: P' = (P * S) / (2 * pi * r_YB^2)
P_prime = (P * S) / (2 * pi * r_YB**2)

print(f"Total power incident on the cell (P'): {P_prime:.4e} W")

# --- Step 4: Convert the answer to microwatts and format ---
P_prime_microwatts = P_prime * 1e6
print("\n--- Final Answer ---")
print(f"The equation for the final power is: P' = (P * S) / (2 * pi * (a - R)^2)")
print(f"P' = ({P:.1e} W * {S} m^2) / (2 * pi * ({a:.4e} m - {R:.4e} m)^2)")
print(f"P' = ({P*S:.1e}) / (2 * pi * ({r_YB:.4e} m)^2)")
print(f"P' = ({P*S:.1e}) / ({2*pi*r_YB**2:.4e})")
print(f"P' = {P_prime:.4e} W")
print(f"The total power incident on the cell is {P_prime_microwatts:.1f} microwatts.")

# Final answer in the required format
final_answer = round(P_prime_microwatts, 1)
# print(f"<<<{final_answer}>>>")
# The above line is commented out to not have it appear in the final output block. The final answer will be hardcoded below.