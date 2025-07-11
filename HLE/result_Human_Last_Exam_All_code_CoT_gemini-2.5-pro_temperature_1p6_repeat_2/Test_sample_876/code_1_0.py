import math

# --- Problem Parameters ---
# Radius of the wires in mm
r_w_mm = 11.25
# Distance from the center of a wire to the center of the cable in mm
m_mm = 17.32
# Inner radius of the common screen in mm
R_mm = 32.32
# Relative permittivity of the insulator
epsilon_r = 4.2

# --- Physical Constants ---
# Permittivity of free space in F/m
epsilon_0 = 8.854e-12

# --- Calculations ---
# The units for radii (mm) cancel out inside the logarithm,
# so we can use them directly for that part of the calculation.
log_argument = (R_mm**2 - m_mm**2) / (r_w_mm * R_mm)

# Denominator of the capacitance formula
denominator = math.log(log_argument)

# Numerator of the capacitance formula
numerator = 2 * math.pi * epsilon_r * epsilon_0

# Capacitance in Farads per meter (F/m)
C_F_per_m = numerator / denominator

# Convert capacitance from F/m to microfarads per kilometer (μF/km)
# 1 F/m = 10^9 μF/km
C_uF_per_km = C_F_per_m * 1e9

# --- Output ---
print("Calculating the capacitance of a three-phase cable with a common screen.")
print("\n--- Formula ---")
print("C (F/m) = (2 * π * ε_r * ε_0) / ln((R² - m²) / (r_w * R))")
print("Result is converted to μF/km by multiplying by 10^9.")

print("\n--- Final Equation with Values ---")
print(f"C (μF/km) = ( (2 * {math.pi:.5f} * {epsilon_r} * {epsilon_0:.5e}) / ln( ( {R_mm}**2 - {m_mm}**2 ) / ( {r_w_mm} * {R_mm} ) ) ) * 10^9")
print(f"C (μF/km) = ( {numerator:.5e} / ln( {log_argument:.5f} ) ) * 10^9")
print(f"C (μF/km) = ( {numerator:.5e} / {denominator:.5f} ) * 10^9")
print(f"C (μF/km) = {C_F_per_m:.5e} * 10^9")

print("\n--- Result ---")
print(f"The capacitance is {C_uF_per_km:.4f} μF/km.")
<<<0.3258>>>