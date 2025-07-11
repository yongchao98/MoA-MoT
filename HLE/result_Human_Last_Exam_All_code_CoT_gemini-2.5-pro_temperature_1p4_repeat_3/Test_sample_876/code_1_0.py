import math

# Step 1: Define constants and parameters in SI units
r_w = 11.25 / 1000  # radius of the wires in meters
m = 17.32 / 1000    # distance from wire center to cable center in meters
R = 32.32 / 1000    # external radius of the cable in meters
epsilon_r = 4.2     # relative permittivity of the insulator (dimensionless)
epsilon_0 = 8.854e-12 # permittivity of free space in F/m

# Step 2: Calculate the terms for the formula
# Numerator of the main formula
numerator = 2 * math.pi * epsilon_0 * epsilon_r

# Argument of the natural logarithm in the denominator
log_argument = (R**3 - m**3) / (3 * m * R * r_w)

# Denominator of the main formula
denominator = math.log(log_argument)

# Step 3: Calculate capacitance in Farads per meter (F/m)
C_F_per_m = numerator / denominator

# Step 4: Convert capacitance to microfarads per kilometer (uF/km)
# 1 F/m = 10^6 uF / (10^-3 km) = 10^9 uF/km
C_uF_per_km = C_F_per_m * 1e9

# Step 5: Print the equation and the final result
print("The formula for capacitance (C) in F/m is:")
print("C = (2 * pi * epsilon_0 * epsilon_r) / ln[(R^3 - m^3) / (3 * m * R * r_w)]\n")

print("Substituting the values in SI units:")
# Using the original values for clarity in the equation string
print(f"C = (2 * {math.pi:.5f} * {epsilon_0} * {epsilon_r}) / ln[({R}^3 - {m}^3) / (3 * {m} * {R} * {r_w})]")

print(f"\nCalculated capacitance per phase:")
print(f"C = {C_uF_per_km:.4f} \u03BCF/km")

# Final answer in the required format
# print(f"<<<{C_uF_per_km:.4f}>>>")