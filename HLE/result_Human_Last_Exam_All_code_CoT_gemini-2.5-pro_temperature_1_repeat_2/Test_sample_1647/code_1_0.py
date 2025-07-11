import math

# This script calculates the time dilation factor 'f' and memory usage 'z'
# according to the problem specification for the Bagua architecture.

# Step 1 & 2: Define physical parameters and calculate the Schwarzschild Radius (Rs).
# Pandora's average orbital radius from Pegasi (m)
pandora_orbital_radius_m = 100_000_000 * 1000
# Pandora's orbital period (s)
pandora_orbital_period_s = 800 * 24 * 3600
# Speed of light (m/s)
speed_of_light_ms = 3e8
# Distance 'd' from the event horizon for the Pioneer probe (m)
pioneer_distance_m = 13 * 1000

# Calculate Rs using the formula derived from Kepler's Third Law.
# Rs = (8 * pi^2 * a^3) / (c^2 * T^2)
rs_numerator = 8 * (math.pi**2) * (pandora_orbital_radius_m**3)
rs_denominator = (speed_of_light_ms**2) * (pandora_orbital_period_s**2)
schwarzschild_radius_rs = rs_numerator / rs_denominator

# Step 3 & 4: Use the binomial approximation for the time dilation factor 'f',
# as the Bagua C environment lacks a sqrt() function.
# f ≈ 1 - Rs / (2 * r), where r = Rs + d
distance_from_center_r = schwarzschild_radius_rs + pioneer_distance_m
time_dilation_factor_f = 1 - (schwarzschild_radius_rs / (2 * distance_from_center_r))

# Round the result to 0.001 as required.
f_rounded = round(time_dilation_factor_f, 3)

# Step 5 & 6: Calculate the memory usage 'z' for the most memory-efficient program.
# The program requires a minimum of 5 variables to store constants and the result:
# 1. a (frac): orbital radius
# 2. T (frac): orbital period
# 3. c (frac): speed of light
# 4. d (unsigned int): probe distance
# 5. result (frac): for storing Rs and then the final factor f
# Each of these data types (frac, int) occupies 8 trits.
num_variables = 5
trits_per_variable = 8
memory_usage_z = num_variables * trits_per_variable

# Step 7: Print the final answer in the format f:z.
# The instruction to "output each number in the final equation" is interpreted
# as printing the two calculated components of the final answer string.
print(f"{f_rounded}:{memory_usage_z}")