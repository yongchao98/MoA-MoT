import math

# Given data in SI units
P = 1e9  # Power of the source in Watts (1 GW)
S = 10.0  # Area of the photovoltaic cell in m^2
T = 12 * 3600  # Orbital period in seconds (12 hours)
A_mirror = 1e6 # Area of the mirror in m^2 (1 km^2)
M = 7.35e22  # Mass of the Moon in kg
R = 1738e3  # Radius of the Moon in meters (1738 km)
G = 6.67e-11  # Gravitational constant in kg^-1 m^3 s^-2

# Step 1: Calculate the semi-major axis 'a' using Kepler's Third Law
# T^2 / a^3 = 4 * pi^2 / (G * M)
# a = (G * M * T^2 / (4 * pi^2))^(1/3)
pi = math.pi
a_cubed = (G * M * T**2) / (4 * pi**2)
a = a_cubed**(1/3)

print(f"Step 1: Calculate orbital properties")
print(f"Gravitational constant G = {G:.2e} kg^-1 m^3 s^-2")
print(f"Lunar mass M = {M:.2e} kg")
print(f"Orbital period T = {T} s")
print(f"Calculated semi-major axis a = {a / 1000:.1f} km")
print("-" * 20)

# Step 2: Determine the distance from the source to the first mirror (d_AX)
# The problem can be solved for any valid eccentricity, so we assume the simplest case: a circular orbit (e=0).
# In a circular orbit, the distance from the center is always 'a'.
# Satellite X is at the zenith of A, so the distance from A to X is the satellite's altitude.
# d_AX = a - R
d_AX = a - R

print(f"Step 2: Analyze geometry and distances")
print(f"Assuming a circular orbit (eccentricity e=0).")
print(f"Lunar radius R = {R / 1000:.1f} km")
print(f"Altitude of the satellite = a - R = {(a - R) / 1000:.1f} km")
print(f"Distance from source A to mirror X, d_AX = {d_AX:.1f} m")
print("-" * 20)

# Step 3: Calculate the final power P' using the virtual source method.
# The derivation shows that the distance from the final virtual source to the cell is d_A''B = 2 * d_AX.
# The power formula is P' = P * S / (4 * pi * (d_A''B)^2)
# P' = P * S / (4 * pi * (2 * d_AX)^2) = P * S / (16 * pi * d_AX^2)
P_prime = (P * S) / (16 * pi * d_AX**2)

print(f"Step 3: Calculate the power on the cell")
print(f"Source power P = {P:.1e} W")
print(f"Cell area S = {S} m^2")
print(f"The final power P' is given by the equation: P' = P * S / (16 * pi * d_AX^2)")
print(f"P' = {P:.1e} * {S} / (16 * {pi:.4f} * {d_AX:.1f}^2)")
print(f"P' = {P * S:.1e} / ({16 * pi * d_AX**2:.2e})")
print(f"Calculated power P' = {P_prime:.4e} W")
print("-" * 20)

# Step 4: Convert the result to microwatts and round to one decimal place.
P_prime_microwatts = P_prime * 1e6

print(f"Step 4: Final Answer")
print(f"The total power incident on the cell is {P_prime_microwatts:.1f} microwatts.")

# Final answer in the required format
final_answer = round(P_prime_microwatts, 1)
# print(f"<<<{final_answer}>>>")