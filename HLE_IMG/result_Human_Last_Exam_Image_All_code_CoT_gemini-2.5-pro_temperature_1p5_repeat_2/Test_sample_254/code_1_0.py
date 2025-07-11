import math

# --- Given Data ---
# Power of the source in Watts
P = 1e9  # 1 GW
# Area of the photovoltaic cell in m^2
S = 10
# Orbital period in seconds
T = 12 * 3600  # 12 hours
# Lunar mass in kg
M = 7.35e22
# Lunar radius in meters
R = 1738e3  # 1738 km
# Gravitational constant in m^3 kg^-1 s^-2
G = 6.67e-11

# --- Step 1: Calculate the semi-major axis (a) of the orbit ---
# Using Kepler's Third Law: T^2 = (4 * pi^2 * a^3) / (G * M)
# We rearrange to solve for a: a = (G * M * T^2 / (4 * pi^2))^(1/3)
pi = math.pi
a_cubed = (G * M * T**2) / (4 * pi**2)
a = a_cubed**(1/3)

# --- Step 2 & 3: Determine the satellite's distance from the Moon's center ---
# We assume a circular orbit (eccentricity e=0) as the simplest case,
# which makes the orbital radius r_X equal to the semi-major axis a.
r_X = a

# --- Step 4 & 5: Calculate the power incident on the cell ---
# The distance from the source A on the surface to satellite X at the zenith.
d_AX = r_X - R

# The source is on the surface, so its power P radiates into a hemisphere (2*pi steradians).
# The intensity at satellite X is I_X = P / (2 * pi * d_AX^2).
# We assume an ideal reflection system where the power density is conserved,
# so the intensity arriving at B is equal to I_X.
# The power on the cell is P_prime = I_X * S.
P_prime_watts = (P * S) / (2 * pi * d_AX**2)

# Convert the power to microwatts (1 W = 1,000,000 uW)
P_prime_microwatts = P_prime_watts * 1e6

# --- Final Output ---
print(f"Calculation Steps:")
print(f"Gravitational constant G = {G} kg^-1 m^3 s^-2")
print(f"Lunar mass M = {M:.2e} kg")
print(f"Orbital period T = {T} s")
print(f"Calculated semi-major axis a = {a:.2f} m")
print(f"Lunar radius R = {R:.2f} m")
print(f"Assuming circular orbit, satellite distance from Moon's center r_X = {r_X:.2f} m")
print(f"Distance from source A to satellite X, d_AX = r_X - R = {d_AX:.2f} m")
print(f"Source power P = {P:.1e} W")
print(f"Cell area S = {S} m^2")
print(f"Final power equation: P' = (P * S) / (2 * pi * d_AX^2)")
print(f"P' = ({P:.1e} * {S}) / (2 * {pi:.4f} * {d_AX:.2f}^2)")
print(f"P' = {P_prime_watts:.3e} W")
print(f"\nFinal Answer:")
print(f"The total power P' incident on the cell is {P_prime_microwatts:.1f} microwatts.")

<<<82.1>>>