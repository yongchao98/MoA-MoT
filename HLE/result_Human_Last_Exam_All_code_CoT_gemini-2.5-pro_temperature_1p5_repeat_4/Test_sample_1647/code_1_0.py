import math

# Plan:
# 1. Calculate the mass (M) of the black hole Pegasi using Kepler's third law
#    M = (4 * pi^2 * a^3) / (G * T^2)
# 2. Calculate the Schwarzschild radius (Rs) of Pegasi.
#    Rs = 2 * G * M / c^2
# 3. Calculate the time dilation factor (f) for the Pioneer probe.
#    f = 1 / sqrt(1 - Rs / r), where r = Rs + d
# 4. Determine the memory usage (z) of a hypothetical C program on the Bagua architecture.
#    This involves summing the memory, in trits, of all required variables.

# --- Step 1 & 2 & 3: Physics Calculations ---

# Constants in SI units
G = 6.67430e-11  # Gravitational constant
C = 299792458    # Speed of light (m/s)

# Pandora's orbital data (convert to SI units)
T_days = 800  # Orbital period in Earth days
T_sec = T_days * 24 * 60 * 60  # Orbital period in seconds

a_km = 100000000  # Average orbital radius in km
a_m = a_km * 1000   # Average orbital radius in meters

# Pioneer probe's distance (convert to SI units)
d_km = 13  # Distance from event horizon in km
d_m = d_km * 1000  # Distance in meters

# Calculate the mass (M) of Pegasi
M_kg = (4 * math.pi**2 * a_m**3) / (G * T_sec**2)

# Calculate the Schwarzschild radius (Rs) of Pegasi
Rs_m = (2 * G * M_kg) / (C**2)

# Calculate the probe's distance (r) from the center of the black hole
r_m = Rs_m + d_m

# Calculate the time dilation factor (f)
# The term inside the square root must be positive.
term = 1 - (Rs_m / r_m)
f = 1 / math.sqrt(term)

# Round f to 3 decimal places
f_rounded = round(f, 3)

# --- Step 4: Memory Usage Calculation (in trits) ---

# Based on the Bagua architecture specification (1 trit = 3 bits)
# - int: 24 bits = 8 trits
# - char: 6 bits = 2 trits
# - frac: 24 bits = 8 trits

# Variables needed in the main function of the C program:
# G, C, PI, T, a, M, Rs, r, term_inside_sqrt, f (10 `frac` types)
# d (1 `int` type)
main_vars_trits = (10 * 8) + (1 * 8)

# A sqrt() function must be implemented. Its local variables also consume memory.
# frac sqrt(frac S) { // S is passed by value, creating a copy
#   frac x;          // Current guess
#   frac temp;        // Temporary for calculation
#   char i;          // Loop counter (0-10 is fine for char)
# }
# S (copy of arg): frac = 8 trits
# x: frac = 8 trits
# temp: frac = 8 trits
# i: char = 2 trits
sqrt_func_vars_trits = 8 + 8 + 8 + 2

# Total memory usage 'z' is the sum of main variables and the sqrt function's stack variables.
z = main_vars_trits + sqrt_func_vars_trits

# --- Final Output ---
# Print the final answer in the required "f:z" format.
# f is the time dilation factor.
# z is the total memory usage in trits.
final_f = f_rounded
final_z = z
print(f"{final_f}:{final_z}")