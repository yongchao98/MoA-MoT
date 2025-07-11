import math

# Given parameters
r_w_mm = 11.25  # Wire radius in mm
m_mm = 17.32    # Distance from wire center to cable center in mm
R_mm = 32.32    # External radius of the cable in mm
eps_r = 4.2     # Relative permittivity of the insulator

# Physical constant
eps_0 = 8.854e-12  # Permittivity of free space in F/m
pi = math.pi

# Convert all lengths from mm to meters for calculation
r_w_m = r_w_mm / 1000
m_m = m_mm / 1000
R_m = R_mm / 1000

# --- Calculation ---

# Calculate the capacitance in Farads per meter (F/m)
numerator = 2 * pi * eps_r * eps_0
log_argument = (R_m**2 - m_m**2) / (R_m * r_w_m)
denominator = math.log(log_argument)
C_F_per_m = numerator / denominator

# Convert the result to microfarads per kilometer (μF/km)
C_uF_per_km = C_F_per_m * 1e9

# --- Output ---

# Print the final equation with all the numbers included as requested
print("The capacitance C is calculated using the following equation (in μF/km):")
print("C = [ (2 * π * ε_r * ε₀) / ln((R² - m²) / (R * r_w)) ] * 10⁹\n")
print("Plugging in the numbers:")
print(
    f"C = [ (2 * {pi:.5f} * {eps_r} * {eps_0}) / "
    f"ln(({R_m}² - {m_m}²) / ({R_m} * {r_w_m})) ] * 10⁹"
)

# Show the evaluated terms for clarity
evaluated_numerator = 2 * pi * eps_r * eps_0
evaluated_log_term = (R_m**2 - m_m**2) / (R_m * r_w_m)
evaluated_log = math.log(evaluated_log_term)
print(f"C = [ {evaluated_numerator:.5e} / ln({evaluated_log_term:.5f}) ] * 10⁹")
print(f"C = [ {evaluated_numerator:.5e} / {evaluated_log:.5f} ] * 10⁹\n")

# Print the final calculated value
print("Final Result:")
print(f"The capacitance of the cable is {C_uF_per_km:.3f} μF/km.")