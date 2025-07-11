import math

# --- User-defined variables from the problem ---
# Radius of the wire in mm
r_w_mm = 11.25
# Distance from wire center to cable center in mm
m_mm = 17.32
# External radius of the cable (screen) in mm
R_mm = 32.32
# Relative permittivity of the insulator
epsilon_r = 4.2

# --- Constants ---
# Permittivity of free space in F/m
epsilon_0 = 8.854e-12

# --- Plan ---
# 1. Calculate absolute permittivity
# 2. Calculate the argument of the natural logarithm
# 3. Calculate capacitance in F/m
# 4. Convert to µF/km and print the result

print("Step 1: Define variables and constants.")
print(f"Wire radius (r_w): {r_w_mm} mm")
print(f"Distance to center (m): {m_mm} mm")
print(f"Screen radius (R): {R_mm} mm")
print(f"Relative permittivity (ε_r): {epsilon_r}")
print(f"Permittivity of free space (ε_0): {epsilon_0:.4e} F/m\n")

# Convert units from mm to m for calculation
r_w = r_w_mm / 1000.0
m = m_mm / 1000.0
R = R_mm / 1000.0

# Calculate absolute permittivity
epsilon = epsilon_r * epsilon_0

# Calculate the argument of the natural logarithm
# Note: The ratio is unitless, so using mm or m gives the same result.
log_arg_num = R**3 - m**3
log_arg_den = 3 * m * R * r_w
log_argument = log_arg_num / log_arg_den

# Calculate the natural logarithm
ln_value = math.log(log_argument)

# Calculate the capacitance in Farads per meter (F/m)
numerator = 2 * math.pi * epsilon
C_per_meter = numerator / ln_value

# Convert capacitance to microfarads per kilometer (µF/km)
# 1 F/m = 1 F/(1e-3 km) = 1e3 F/km = 1e3 * 1e6 µF/km = 1e9 µF/km
C_per_km = C_per_meter * 1e9

print("Step 2: Show the calculation using the formula.")
print("The formula is: C = (2 * π * ε_r * ε_0) / ln((R³ - m³) / (3 * m * R * r_w))\n")

print("--- Final Equation with numerical values ---")
# The final equation is constructed by substituting all the numerical values.
# To keep it readable, we show the geometric parts in mm, as the ratio is dimensionless.
final_equation_str = (
    f"C (µF/km) = (2 * {math.pi:.5f} * {epsilon_r} * {epsilon_0:.4e} F/m) / "
    f"ln(({R_mm}³ - {m_mm}³) / (3 * {m_mm} * {R_mm} * {r_w_mm})) * 10⁹"
)
print(final_equation_str)
print("-" * 35)

# Calculate intermediate values for clarity
num_calc = R_mm**3 - m_mm**3
den_calc = 3 * m_mm * R_mm * r_w_mm
print(f"Numerator of ln: {R_mm}³ - {m_mm}³ = {num_calc:.2f}")
print(f"Denominator of ln: 3 * {m_mm} * {R_mm} * {r_w_mm} = {den_calc:.2f}")
print(f"ln argument: {num_calc:.2f} / {den_calc:.2f} = {log_argument:.5f}")
print(f"ln({log_argument:.5f}) = {ln_value:.5f}\n")
print(f"Capacitance in µF/km = {C_per_km:.4f}")

<<<0.5617>>>