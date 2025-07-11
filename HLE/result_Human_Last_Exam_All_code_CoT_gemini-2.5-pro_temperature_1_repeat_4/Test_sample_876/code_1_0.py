import math

# Step 1: Define the given parameters and constants
# Given values in millimeters
r_w_mm = 11.25  # Wire radius in mm
m_mm = 17.32    # Distance from wire center to cable center in mm
R_mm = 32.32    # External radius of the cable in mm

# Given relative permittivity
epsilon_r = 4.2

# Physical constant: permittivity of free space in F/m
epsilon_0 = 8.854e-12

# Step 2: Convert dimensions from mm to meters for calculation
r_w = r_w_mm / 1000.0
m = m_mm / 1000.0
R = R_mm / 1000.0

# Step 3: Calculate the capacitance
# Absolute permittivity of the insulator
epsilon = epsilon_r * epsilon_0

# Calculate the numerator and denominator of the argument of the logarithm
log_arg_numerator = R**3 - m**3
log_arg_denominator = 3 * m * R * r_w

# Calculate the argument
log_argument = log_arg_numerator / log_arg_denominator

# Calculate capacitance in Farads per meter (F/m)
C_per_meter = (2 * math.pi * epsilon) / math.log(log_argument)

# Step 4: Convert the result to microfarads per kilometer (μF/km)
C_per_km_uF = C_per_meter * 1e9

# Step 5: Print the equation with values and the final result
print("The formula for capacitance per phase is C = (2 * pi * epsilon) / ln( (R^3 - m^3) / (3 * m * R * r_w) )")
print("\nSubstituting the values (in SI units):")

equation_str = (
    f"C = (2 * {math.pi:.5f} * {epsilon:.4e}) / "
    f"ln( ({R}**3 - {m}**3) / (3 * {m} * {R} * {r_w}) )"
)
print(equation_str)

print(f"\nCalculated value of (R^3 - m^3) = {log_arg_numerator:.4e}")
print(f"Calculated value of (3 * m * R * r_w) = {log_arg_denominator:.4e}")
print(f"Calculated value of ln(...) = {math.log(log_argument):.4f}")

print(f"\nThe capacitance is {C_per_meter:.4e} F/m.")
print(f"\nFinal Answer: The capacitance of the cable is {C_per_km_uF:.4f} uF/km.")
print(f"<<<{C_per_km_uF:.4f}>>>")