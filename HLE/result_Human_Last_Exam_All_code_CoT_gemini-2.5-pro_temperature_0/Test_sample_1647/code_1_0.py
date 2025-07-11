import math

# This script calculates the gravitational time dilation factor (f) and
# the memory usage (z) for a C program on the Bagua architecture.

# Step 1: Calculate the Schwarzschild radius (Rs) of the black hole Pegasi.
# Constants
G = 6.674e-11  # Gravitational constant in N(m/kg)^2
C = 3e8        # Speed of light in m/s

# Pandora's orbital data (convert to SI units)
a_pandora_km = 100_000_000
a_pandora_m = a_pandora_km * 1000  # Orbital radius in meters

t_pandora_days = 800
t_pandora_s = t_pandora_days * 24 * 3600  # Orbital period in seconds

# Mass of Pegasi (M) from Kepler's Third Law: M = (4 * pi^2 * a^3) / (G * T^2)
M_pegasi = (4 * math.pi**2 * a_pandora_m**3) / (G * t_pandora_s**2)

# Schwarzschild Radius (Rs): Rs = 2 * G * M / c^2
Rs_meters = (2 * G * M_pegasi) / (C**2)

# Step 2: Calculate the time dilation factor (f).
# Given distance from the event horizon
d_km = 13
d_meters = d_km * 1000

# Time dilation factor: f = sqrt(d / (Rs + d))
factor = math.sqrt(d_meters / (Rs_meters + d_meters))
f_rounded = round(factor, 3)

# Step 3: Determine the memory usage (z) in trits.
# The most memory-efficient C program would define Rs as a constant and use
# variables for the input 'd' and the fractional result.
# unsigned int d;  // 8 trits
# frac result;     // 8 trits
z_trits = 8 + 8

# Step 4: Output the final equation and the answer.
# As required, here is the final equation with each number:
print("Final Equation:")
print(f"f = sqrt(d / (Rs + d))")
print(f"f = sqrt({d_meters} / ({Rs_meters:.0f} + {d_meters}))")
print("")

# Print the final answer in the specified f:z format.
print("Answer:")
print(f"{f_rounded}:{z_trits}")