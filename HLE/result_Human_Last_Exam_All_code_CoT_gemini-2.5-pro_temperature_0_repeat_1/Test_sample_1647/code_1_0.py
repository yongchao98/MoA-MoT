import math

# This script calculates the gravitational time dilation factor 'f' for the Pioneer probe
# and determines the memory usage 'z' of an efficient C program for the Bagua architecture.

# --- Step 1 & 2: Calculate Schwarzschild Radius (Rs) ---
# We use a combined formula derived from Kepler's Third Law and the Rs formula:
# Rs = (8 * pi^2 * a^3) / (c^2 * T^2)

# Constants from the problem description, converted to SI units.
# Pandora's average orbital radius in meters.
a = 100000000 * 1000
# Pandora's orbital period in seconds.
T = 800 * 24 * 3600
# Speed of light in m/s.
c = 299792458

# Perform the calculation for Rs.
Rs_numerator = 8 * (math.pi**2) * (a**3)
Rs_denominator = (c**2) * (T**2)
Rs = Rs_numerator / Rs_denominator

# --- Step 3: Calculate Time Dilation Factor (f) ---
# The simplified formula for time dilation is f = sqrt(1 + Rs / d).

# Pioneer's distance from the event horizon in meters.
d = 13 * 1000

# Calculate the value of f.
f_value = math.sqrt(1 + Rs / d)

# Round f to 3 decimal places (0.001).
f_rounded = round(f_value, 3)

# --- Step 4: Determine Memory Usage (z) ---
# A memory-efficient C program for Bagua would declare only essential variables.
# The Bagua architecture specifies that 'int' and 'frac' types are both 8 trits.
# 1. An 'int' to store the distance 'd' (13000). Size: 8 trits.
# 2. A 'frac' to store the calculated Schwarzschild radius 'Rs'. Size: 8 trits.
# 3. A 'frac' to store the final factor 'f'. Size: 8 trits.
num_variables = 3
trits_per_variable = 8
z_memory_usage = num_variables * trits_per_variable

# --- Step 5: Output the final result ---
# The problem asks for the answer in the format f:z.
# The final equation is f:z, and the numbers are the calculated f_rounded and z_memory_usage.
final_f = f_rounded
final_z = z_memory_usage

print(f"{final_f}:{final_z}")