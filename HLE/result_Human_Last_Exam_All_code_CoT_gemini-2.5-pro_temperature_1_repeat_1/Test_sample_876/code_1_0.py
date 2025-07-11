import math

# Description of the task and formula used
print("This script calculates the capacitance of a three-phase cable with a common screen.")
print("The formula used is: C = (2 * pi * ε_r * ε_0) / ln((R^3 - m^3) / (3 * r_w * m * R))")
print("-" * 30)

# Given parameters
r_w = 11.25  # wire radius in mm
m = 17.32    # distance from wire center to cable center in mm
R = 32.32    # external radius of the cable in mm
eps_r = 4.2  # relative permittivity of the insulator

# Physical constant
eps_0 = 8.854e-12  # Permittivity of free space in F/m

# --- Calculation ---

# The argument of the logarithm is dimensionless, so we can use mm directly.
log_arg_numerator = R**3 - m**3
log_arg_denominator = 3 * r_w * m * R
log_argument = log_arg_numerator / log_arg_denominator

# Calculate the capacitance in Farads per meter (F/m)
capacitance_F_per_m = (2 * math.pi * eps_r * eps_0) / math.log(log_argument)

# Convert the capacitance to microfarads per kilometer (µF/km)
# 1 F/m = 1e9 µF/km
capacitance_muF_per_km = capacitance_F_per_m * 1e9

# --- Output the result ---

print("Equation with the given values:")
# Printing each number in the equation as requested
equation_str = (
    f"C = (2 * {math.pi:.5f} * {eps_r} * {eps_0}) / "
    f"ln(({R}**3 - {m}**3) / (3 * {r_w} * {m} * {R}))"
)
print(equation_str)
print("-" * 30)
print(f"The calculated capacitance is: {capacitance_muF_per_km:.3f} µF/km")
print(f"<<<{capacitance_muF_per_km:.3f}>>>")