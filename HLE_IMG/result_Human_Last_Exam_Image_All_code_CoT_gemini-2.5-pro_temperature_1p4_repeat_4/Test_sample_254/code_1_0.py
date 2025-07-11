import math

# Given data
P = 1e9  # Luminosity of the source in Watts (1 GW)
S_cell = 10  # Area of the photovoltaic cell in m^2
A_mirror = 1e6 # Area of the satellite mirror in m^2 (1 km^2)

M = 7.35e22  # Mass of the Moon in kg
R = 1738e3   # Radius of the Moon in m (1738 km)
G = 6.67e-11 # Gravitational constant in m^3 kg^-1 s^-2
T = 12 * 3600 # Orbital period in seconds (12 hours)

# Step 1: Calculate the semi-major axis 'a' using Kepler's Third Law
# T^2 / a^3 = (4 * pi^2) / (G * M)
# a = (G * M * T^2 / (4 * pi^2))^(1/3)
four_pi_squared = 4 * math.pi**2
a_cubed = (G * M * T**2) / four_pi_squared
a = a_cubed**(1/3)

# Step 2: Calculate the altitude of satellite X
# We assume the relevant distance r_X is the semi-major axis a.
# The distance from the source A on the surface to satellite X is its altitude.
d_AX = a - R

# Step 3: Calculate the power P' on the cell
# The simplified formula assumes the geometric factors cancel out to 1.
# P' = (P * S_cell) / (4 * pi * d_AX^2)
P_prime = (P * S_cell) / (4 * math.pi * d_AX**2)

# Step 4: Convert the power to microwatts
P_prime_microwatts = P_prime * 1e6

# --- Output the results with explanations ---
print("--- Calculation Steps ---")
print(f"Given Luminosity of Source (P): {P:.2e} W")
print(f"Given Area of Photovoltaic Cell (S): {S_cell} m^2")
print(f"Given Mass of the Moon (M): {M:.2e} kg")
print(f"Given Radius of the Moon (R): {R/1e3:.0f} km")
print(f"Given Orbital Period (T): {T/3600:.0f} hours")
print("\nStep 1: Calculate the semi-major axis (a)")
print(f"a = (G * M * T^2 / (4 * pi^2))^(1/3)")
print(f"a = ({G:.2e} * {M:.2e} * {T**2:.2e} / (4 * {math.pi**2:.2f}))^(1/3)")
print(f"a = {a/1e3:.1f} km")
print("\nStep 2: Calculate the altitude of satellite X (d_AX)")
print(f"d_AX = a - R")
print(f"d_AX = {a/1e3:.1f} km - {R/1e3:.1f} km = {d_AX/1e3:.1f} km")
print("\nStep 3: Calculate the power incident on the cell (P')")
print("Using the simplified power transfer formula: P' = (P * S_cell) / (4 * pi * d_AX^2)")
print(f"P' = ({P:.2e} W * {S_cell} m^2) / (4 * pi * ({d_AX:.2e} m)^2)")
print(f"P' = {P_prime:.4e} W")
print("\nStep 4: Convert power to microwatts")
print(f"P' (uW) = {P_prime:.4e} W * 1,000,000 uW/W")
print(f"P' = {P_prime_microwatts:.1f} uW")
print("\n--- Final Answer ---")
print(f"The total power P' incident on the cell is {P_prime_microwatts:.1f} microwatts.")
