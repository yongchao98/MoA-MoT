import math

# --- Given parameters ---
# Radius of the wire in mm
r_w_mm = 11.25
# Distance from wire center to cable center in mm
m_mm = 17.32
# External radius of the cable in mm
R_mm = 32.32
# Relative permittivity of the insulator
epsilon_r = 4.2
# Permittivity of free space in F/m
epsilon_0 = 8.854e-12

# --- Calculation ---
# The units (mm) cancel out inside the logarithm argument
# Calculate R^3 and m^3
R_cubed = R_mm**3
m_cubed = m_mm**3

# Calculate the argument of the natural logarithm
log_argument = (R_mm / r_w_mm) * (R_cubed - m_cubed) / (R_cubed + m_cubed)

# Calculate the denominator of the main formula
denominator = math.log(log_argument)

# Calculate the numerator of the main formula
numerator = 2 * math.pi * epsilon_r * epsilon_0

# Calculate the capacitance in Farads per meter (F/m)
capacitance_F_per_m = numerator / denominator

# Convert from F/m to microfarads per kilometer (uF/km)
# Conversion factor: 1 F/m = 10^9 uF/km
conversion_factor = 1e9
capacitance_uF_per_km = capacitance_F_per_m * conversion_factor

# --- Output ---
print("This script calculates the capacitance of a three-phase cable with a common screen.")
print("\nFormula used for capacitance to neutral (C):")
print("C = (2 * pi * ε_r * ε_0) / ln[ (R/r_w) * (R^3 - m^3) / (R^3 + m^3) ]\n")

print("Given values:")
print(f"Wire radius (r_w) = {r_w_mm} mm")
print(f"Distance (m) = {m_mm} mm")
print(f"External radius (R) = {R_mm} mm")
print(f"Relative permittivity (ε_r) = {epsilon_r}")
print(f"Permittivity of free space (ε_0) = {epsilon_0} F/m\n")

print("Final equation with values (for result in μF/km):")
# We use the raw numbers in the equation string for clarity
equation_str = (
    f"C (μF/km) = (2 * π * {epsilon_r} * {epsilon_0:.4e}) / "
    f"ln[({R_mm}/{r_w_mm}) * ({R_mm}³ - {m_mm}³)/({R_mm}³ + {m_mm}³)] * {conversion_factor:.0e}"
)
print(equation_str)

print(f"\nCalculated Capacitance = {capacitance_uF_per_km:.4f} μF/km")
