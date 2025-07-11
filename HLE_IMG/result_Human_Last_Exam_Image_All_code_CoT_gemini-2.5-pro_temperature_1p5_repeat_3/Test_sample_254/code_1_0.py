import math

# Given constants
P = 1e9  # W
S = 10.0  # m^2
G = 6.67e-11  # m^3 kg^-1 s^-2
M = 7.35e22  # kg
R = 1738e3  # m
T = 12 * 3600  # s

# 1. Calculate the semi-major axis 'a'
# T^2 / a^3 = 4 * pi^2 / (G * M)
# a^3 = G * M * T^2 / (4 * pi^2)
mu = G * M
a_cubed = mu * T**2 / (4 * math.pi**2)
a = a_cubed**(1/3)

# 2. Calculate the power P' incident on the cell
# The derived formula for power is P' = P * S / (2 * pi * h_a^2)
# With the deduced constraint that the apogee altitude h_a = a
h_a = a
P_prime = (P * S) / (2 * math.pi * h_a**2)

# 3. Convert to microwatts and print
P_prime_microwatts = P_prime * 1e6

print(f"Gravitational parameter (G*M) = {mu:.4e} m^3/s^2")
print(f"Orbital period T = {T} s")
print(f"Semi-major axis a = {a:.4f} m")
# This is our key assumption for solvability
print(f"Apogee altitude h_a is assumed to be equal to a, h_a = {h_a:.4f} m")
print(f"\nFinal power equation:")
print(f"P' = P * S / (2 * pi * a^2)")
print(f"P' = {P:.0e} W * {S} m^2 / (2 * pi * ({a:.4f} m)^2)")
print(f"P' = {P_prime:.4e} W")
print(f"P' = {P_prime_microwatts:.1f} microwatts")
