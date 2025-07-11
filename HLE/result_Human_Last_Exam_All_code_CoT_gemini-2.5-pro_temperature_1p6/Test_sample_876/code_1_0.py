import math

# Description:
# This script calculates the capacitance per phase of a three-phase cable
# with a common screen.

# --- Given values and constants ---

# Wire radius in millimeters and meters
r_w_mm = 11.25
r_w_m = r_w_mm / 1000

# Distance from wire center to cable center in millimeters and meters
m_mm = 17.32
m_m = m_mm / 1000

# External radius of the cable (internal radius of the screen) in millimeters and meters
R_mm = 32.32
R_m = R_mm / 1000

# Relative permittivity of the insulator (dimensionless)
epsilon_r = 4.2

# Permittivity of free space in Farads per meter (F/m)
epsilon_0 = 8.854e-12

# --- Calculation ---

# Calculate the numerator of the capacitance formula
numerator = 2 * math.pi * epsilon_0 * epsilon_r

# Calculate the argument of the natural logarithm in the denominator
log_argument = (R_m**2 - m_m**2) / (R_m * r_w_m)

# Calculate the denominator
denominator = math.log(log_argument)

# Calculate capacitance in Farads per meter (F/m)
C_F_per_m = numerator / denominator

# Convert capacitance to microfarads per kilometer (μF/km)
# 1 F/m = 10^9 μF/km
C_uF_per_km = C_F_per_m * 1e9

# --- Output the result ---

print("The formula for capacitance (C) is: C = (2 * pi * epsilon_0 * epsilon_r) / ln((R^2 - m^2) / (R * r_w))")
print("\nSubstituting the given values into the formula (in SI units):")
# The final code prints the full equation with each number.
# Using f-string formatting to embed the numeric values directly into the equation string.
print(f"C = (2 * {math.pi:.6f} * {epsilon_0} * {epsilon_r}) / ln(({R_m}**2 - {m_m}**2) / ({R_m} * {r_w_m}))")
print(f"\nThe calculated capacitance is: {C_uF_per_km:.4f} μF/km")

<<<0.3260>>>