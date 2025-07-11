import math

# This script calculates the gravitational time dilation factor 'f' for the Pioneer probe
# and the memory usage 'z' for an optimized program on the Bagua architecture.

# Part 1: Calculate the time dilation factor 'f'

# Define physical constants in SI units
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
C = 299792458    # Speed of light (m/s)

# Pandora's orbital parameters converted to SI units
R_ORB_KM = 100_000_000
r_orb_m = R_ORB_KM * 1000  # Orbital radius in meters
T_DAYS = 800
t_s = T_DAYS * 24 * 60 * 60  # Orbital period in seconds

# Pioneer's distance from the event horizon
d_km = 13

# Calculate the mass of the black hole Pegasi (M) using Kepler's Third Law
# T^2 = (4 * pi^2 * R^3) / (G * M)  =>  M = (4 * pi^2 * R^3) / (G * T^2)
pegasi_mass = (4 * math.pi**2 * r_orb_m**3) / (G * t_s**2)

# Calculate the Schwarzschild Radius (Rs) for Pegasi
# Rs = 2 * G * M / c^2
rs_m = (2 * G * pegasi_mass) / (C**2)
rs_km = rs_m / 1000

# Calculate the gravitational time dilation factor (f)
# f = sqrt(d / (Rs + d))
# All units in the formula must be consistent (km)
time_dilation_factor = math.sqrt(d_km / (rs_km + d_km))

# Round the factor to 3 decimal places as required
f_rounded = round(time_dilation_factor, 3)


# Part 2: Calculate the memory usage 'z'

# The most memory-efficient C program for Bagua would only declare variables
# for the essential constants and the final result.
# The value d=13 can be used as a literal in the calculation to save memory.
# - One 'frac' variable is needed to store the Schwarzschild Radius (Rs).
# - One 'frac' variable is needed for the iterative calculation of the square root (f).
#
# From the Bagua specification, the size of a 'frac' is:
# size(signed char) + size(unsigned wchar) + size(signed char)
# 2 trits + 4 trits + 2 trits = 8 trits
sizeof_frac_trits = 8
num_variables = 2
z_memory_usage = sizeof_frac_trits * num_variables


# Part 3: Print the final answer in the format 'f:z'
# The final equation is the requested answer format 'f:z'.
# The numbers in this equation are the calculated values of 'f' and 'z'.
final_f = f_rounded
final_z = z_memory_usage

print(f"{final_f}:{final_z}")