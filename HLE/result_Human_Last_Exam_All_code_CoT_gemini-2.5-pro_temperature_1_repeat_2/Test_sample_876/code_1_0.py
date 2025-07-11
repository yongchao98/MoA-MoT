import math

# --- 1. Define variables and constants ---

# Given values
r_w_mm = 11.25  # Wire radius in mm
m_mm = 17.32    # Distance from wire center to cable center in mm
R_mm = 32.32    # External radius of the cable (inner radius of screen) in mm
epsilon_r = 4.2 # Relative permittivity of the insulator

# Convert dimensions from mm to m for calculation
r_w = r_w_mm / 1000
m = m_mm / 1000
R = R_mm / 1000

# Physical constant
epsilon_0 = 8.854e-12  # Permittivity of free space in F/m

# --- 2. Calculate capacitance using the formula ---

# Numerator of the main formula
numerator_C = 2 * math.pi * epsilon_0 * epsilon_r

# Argument of the natural logarithm
log_argument_num = (R**3 - m**3)
log_argument_den = (3 * m * R * r_w)
log_argument = log_argument_num / log_argument_den

# Denominator of the main formula
denominator_C = math.log(log_argument)

# Capacitance in Farads per meter (F/m)
C_F_per_m = numerator_C / denominator_C

# --- 3. Convert to microfarads per kilometer (uF/km) ---
# 1 F/m = 1e6 uF / 1e-3 km = 1e9 uF/km
C_uF_per_km = C_F_per_m * 1e9

# --- 4. Print the output ---

print("The formula for capacitance per phase for a three-phase cable with a common screen is:")
print("C = (2 * pi * epsilon_0 * epsilon_r) / ln((R^3 - m^3) / (3 * m * R * r_w))")
print("\nPlugging in the values (in SI units):")

# Using an f-string to format the output equation clearly
equation = (
    f"C (F/m) = (2 * {math.pi:.5f} * {epsilon_0} * {epsilon_r}) / "
    f"ln((({R})^3 - ({m})^3) / (3 * {m} * {R} * {r_w}))"
)
print(equation)

print(f"\nCalculated Capacitance = {C_uF_per_km:.4f} \u03BCF/km")

# Final answer in the required format
final_answer = f"<<<{C_uF_per_km:.4f}>>>"
# This print statement is for the platform, the user sees the one above.
# print(final_answer)