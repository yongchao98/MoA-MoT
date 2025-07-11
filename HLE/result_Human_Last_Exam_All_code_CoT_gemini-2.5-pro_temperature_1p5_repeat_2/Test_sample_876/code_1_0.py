import math

# --- Task: Calculate the capacitance of a three-phase cable ---

# --- Given parameters ---
r_w_mm = 11.25  # Radius of the wire in mm
m_mm = 17.32    # Distance from wire center to cable center in mm
R_mm = 32.32    # External radius of the cable (screen radius) in mm
eps_r = 4.2     # Relative permittivity of the insulator

# --- Constants ---
# Permittivity of free space in Farads per meter (F/m)
eps_0 = 8.854e-12

# The calculation requires the capacitance in microfarads per kilometer (μF/km).
# The conversion factor from F/m to μF/km is 10^9.
# (1 F/m) * (10^6 μF / 1 F) * (1000 m / 1 km) = 10^9 μF/km
CONVERSION_FACTOR = 1e9

# --- Calculation ---
# The formula for capacitance C in μF/km is:
# C = (2 * π * ε_r * ε₀ * 10⁹) / ln((R³ - m³) / (3 * m * R * r_w))
# Note: The units for R, m, and r_w inside the logarithm cancel out,
# so we can use them directly in mm.

# Calculate the numerator of the formula, including the conversion factor
numerator = 2 * math.pi * eps_r * eps_0 * CONVERSION_FACTOR

# Calculate the argument of the natural logarithm (ln)
log_argument = (R_mm**3 - m_mm**3) / (3 * m_mm * R_mm * r_w_mm)

# Calculate the denominator of the formula
denominator = math.log(log_argument)

# Calculate the final capacitance
capacitance = numerator / denominator

# --- Output the results ---
# As requested, printing the final equation with the numbers plugged in.
print("The formula for capacitance (C) in μF/km is:")
print("C = (2 * π * ε_r * ε₀ * 10⁹) / ln((R³ - m³) / (3 * m * R * r_w))")
print("\nSubstituting the given values into the equation:")
print(f"C = (2 * 3.14159... * {eps_r} * {eps_0} * 10^9) / ln(({R_mm}³ - {m_mm}³) / (3 * {m_mm} * {R_mm} * {r_w_mm}))")

# Print the evaluated numerator and denominator for clarity
print("\nEvaluated terms:")
print(f"Numerator = {numerator:.4f}")
print(f"Denominator = ln({log_argument:.4f}) = {denominator:.4f}")

# Print the final result
print("\nFinal Calculation:")
print(f"C = {numerator:.4f} / {denominator:.4f}")
print(f"The capacitance of the cable is {capacitance:.4f} μF/km.")
<<<0.5634>>>