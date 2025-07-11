import math

# Step 1 & 2: Calculate the Schwarzschild Radius (Rs) of Pegasi
# Constants
G = 6.674e-11  # Gravitational constant
c = 299792458  # Speed of light in m/s
pi = math.pi

# Pandora's orbital data
a_km = 100_000_000  # Orbital radius in km
T_days = 800        # Orbital period in Earth days

# Convert to SI units (meters and seconds)
a_m = a_km * 1000
T_s = T_days * 24 * 60 * 60

# Calculate the mass of the black hole Pegasi (M) using Kepler's Third Law
# M = (4 * pi^2 * a^3) / (G * T^2)
M_pegasi = (4 * pi**2 * a_m**3) / (G * T_s**2)

# Calculate the Schwarzschild Radius (Rs) in meters
Rs_m = (2 * G * M_pegasi) / (c**2)

# Convert Rs to kilometers for the time dilation formula
Rs_km = Rs_m / 1000

# Step 3: Calculate the Time Dilation Factor (f)
# Given distance from the event horizon (d)
d_km = 13  # in km

# Formula: f = sqrt(1 + Rs / d)
f_val = math.sqrt(1 + (Rs_km / d_km))

# Round f to 3 decimal places (0.001)
f_rounded = round(f_val, 3)

# Step 4: Calculate the memory usage (z) for the Bagua C program
# Based on the most memory-efficient C program design:
# A sqrt_frac function taking 'frac x' and declaring 'frac guess', 'int i'.
# The main function declares 'frac f'.
# Total variables: 3 of type 'frac', 1 of type 'int'.
frac_size_trits = 8
int_size_trits = 8
num_frac_vars = 3  # f, x, guess
num_int_vars = 1   # i

z_memory_usage = (num_frac_vars * frac_size_trits) + (num_int_vars * int_size_trits)

# Step 5: Output the final answer
# As requested, output the numbers in the final equation first.
print(f"Calculating the time dilation factor 'f' for the given distance d = {d_km} km:")
print(f"Using the formula: f = sqrt(1 + Rs/d)")
print(f"Plugging in the values: f = sqrt(1 + {Rs_km:.5f}/{d_km})")
print(f"Result for f (rounded to 0.001): {f_rounded}")
print("-" * 20)
print(f"Calculating memory usage 'z' for the C program variables:")
print(f"Variables: 3 'frac' (size {frac_size_trits} trits each), 1 'int' (size {int_size_trits} trits)")
print(f"z = (3 * {frac_size_trits}) + (1 * {int_size_trits}) = {z_memory_usage} trits")
print("-" * 20)
print("Final Answer (f:z):")
print(f"{f_rounded}:{z_memory_usage}")