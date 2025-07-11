import math

# Define the given values and constants
# Radius of the wire in mm
r_w = 11.25
# Distance from the center of a wire to the center of the cable in mm
m = 17.32
# External radius of the cable in mm
R = 32.32
# Permittivity of the insulator (relative permittivity)
epsilon = 4.2
# Permittivity of free space in F/m
epsilon_0 = 8.854e-12

# The argument of the natural logarithm is dimensionless,
# so we can use the values in mm directly inside the log function.
log_numerator = R**3 - m**3
log_denominator = 3 * m * r_w * R
log_argument = log_numerator / log_denominator

# Calculate the capacitance in Farads per meter (F/m)
# The formula is C = (2 * pi * epsilon_0 * epsilon) / ln(log_argument)
numerator = 2 * math.pi * epsilon_0 * epsilon
denominator = math.log(log_argument)
capacitance_per_meter = numerator / denominator

# Convert capacitance from F/m to µF/km
# 1 F/m = 1000 F/km
# 1 F/km = 1e6 µF/km
# So, conversion factor is 1000 * 1e6 = 1e9
capacitance_uF_per_km = capacitance_per_meter * 1e9

# Print the final equation with all the numbers
# Note: The conversion factor 1e9 is used to get the result in µF/km
print("The capacitance C in µF/km is calculated using the formula:")
print(f"C = (2 * pi * {epsilon_0} * {epsilon}) / ln(({R}^3 - {m}^3) / (3 * {m} * {r_w} * {R})) * 1e9")
print("\nFinal Calculation:")
print(f"C = {capacitance_uF_per_km:.4f} µF/km")