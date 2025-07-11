import math

# Step 1: Define the parameters and constants
# Given parameters from the user
r_w_mm = 11.25  # radius of the wires in mm
m_mm = 17.32    # distance from wire center to cable center in mm
R_mm = 32.32    # external radius of the cable in mm
epsilon_r = 4.2 # relative permittivity of the insulator

# Physical constant
epsilon_0 = 8.854e-12  # Permittivity of free space in F/m

# Step 2: Convert units from mm to meters
r_w = r_w_mm / 1000.0
m = m_mm / 1000.0
R = R_mm / 1000.0

# Step 3: Calculate the permittivity of the insulator
epsilon = epsilon_r * epsilon_0

# Step 4: Calculate the capacitance using the formula
# The formula for capacitance per phase (in F/m) is:
# C = (2 * pi * epsilon) / ln((R^2 - m^2) / (R * r_w))

# Calculate the numerator and the terms in the denominator
numerator = 2 * math.pi * epsilon
log_argument = (R**2 - m**2) / (R * r_w)
denominator = math.log(log_argument)

# Calculate capacitance in F/m
C_per_meter = numerator / denominator

# Step 5: Convert the result to microfarads per kilometer (uF/km)
# 1 F/m = 10^9 uF/km
C_per_km_uF = C_per_meter * 1e9

# Step 6: Print the results in a clear and detailed format
print("--- Calculation of Three-Phase Cable Capacitance ---")
print("\nGiven Parameters:")
print(f"Wire radius (r_w): {r_w_mm} mm = {r_w} m")
print(f"Wire-to-center distance (m): {m_mm} mm = {m} m")
print(f"External radius (R): {R_mm} mm = {R} m")
print(f"Relative permittivity of insulator (ε_r): {epsilon_r}")

print("\nFormula Used:")
print("C = (2 * π * ε) / ln((R² - m²) / (R * r_w))")

print("\nEquation with numerical values plugged in:")
print(f"C (F/m) = (2 * {math.pi:.5f} * {epsilon:.5e}) / ln(({R}**2 - {m}**2) / ({R} * {r_w}))")
print(f"C (F/m) = ({numerator:.5e}) / ln({log_argument:.5f})")
print(f"C (F/m) = ({numerator:.5e}) / {denominator:.5f}")
print(f"C = {C_per_meter:.5e} F/m")

print("\n--- Final Result ---")
print(f"The capacitance of the cable is {C_per_km_uF:.4f} μF per km.")

# For final answer extraction
final_answer = round(C_per_km_uF, 4)
# print(f'<<<{final_answer}>>>')