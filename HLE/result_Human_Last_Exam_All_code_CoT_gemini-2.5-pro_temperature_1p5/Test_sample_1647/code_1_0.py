import math

# Step 1 & 2: Define constants and simplify calculations
# Input data from the problem description
pandora_orbital_radius_km = 100_000_000
pandora_orbital_period_days = 800
pioneer_distance_from_horizon_km = 13

# Constants
SECS_PER_DAY = 86400
METERS_PER_KM = 1000
C = 299792458  # Speed of light in m/s

# Convert inputs to SI units (meters and seconds)
a_m = pandora_orbital_radius_km * METERS_PER_KM
T_sec = pandora_orbital_period_days * SECS_PER_DAY
d_m = pioneer_distance_from_horizon_km * METERS_PER_KM

# Pre-calculate squares
PI2 = math.pi**2
C2 = C**2
T_sec2 = T_sec**2

# Simplified formula for Schwarzschild Radius (Rs) where G cancels out:
# Rs = (8 * pi^2 * a^3) / (c^2 * T^2)
Rs = (8 * PI2 * (a_m**3)) / (C2 * T_sec2)

# Simplified formula for time dilation factor (f), using approximation:
# f ≈ 1 - Rs / (2 * (Rs + d))
# This is valid because Rs is much smaller than d.
f = 1 - (Rs / (2 * (Rs + d_m)))

# Round f to 3 decimal places
f_rounded = round(f, 3)

# Step 3 & 4: Determine memory usage for the most efficient Bagua C program
# The simplified calculation requires the following variables in the C program:
# 1. a_m: Pandora's orbital radius in meters (as frac)
# 2. T_sec: Pandora's orbital period in seconds (as frac)
# 3. d_m: Pioneer's distance in meters (as frac)
# 4. PI2: The constant pi squared (as frac)
# 5. C2: The constant c squared (as frac)
# 6. Rs: The calculated Schwarzschild radius (as frac)
# 7. f: The final time dilation factor (as frac)
#
# Total number of variables = 7
# Each 'frac' variable occupies 8 trits.
num_variables = 7
trits_per_variable = 8
z_memory_usage = num_variables * trits_per_variable

# Step 5: Output the result in the specified format
# The final equation is the f:z representation.
# The code prints each part of that "equation": the number f, the colon, and the number z.
print(f"{f_rounded}:{z_memory_usage}")