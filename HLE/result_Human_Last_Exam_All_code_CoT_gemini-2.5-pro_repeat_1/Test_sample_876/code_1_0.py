import math

# --- Task: Calculate the capacitance of a three-phase cable with a common screen ---

# --- Given values ---
r_w_mm = 11.25  # radius of the wire in mm
m_mm = 17.32    # distance from wire center to cable center in mm
R_mm = 32.32    # external radius of the cable (screen) in mm
epsilon_r = 4.2 # relative permittivity of the insulator

# --- Constants ---
epsilon_0 = 8.854e-12 # permittivity of free space in F/m

# --- Plan ---
# 1. Convert all length units to SI units (meters).
# 2. Calculate the term inside the natural logarithm.
# 3. Calculate the full capacitance formula in Farads per meter (F/m).
# 4. Convert the final result to microfarads per kilometer (uF/km).
# 5. Print the final equation with all numbers and the result.

# --- Step 1: Convert units to SI (meters) ---
r_w = r_w_mm / 1000
m = m_mm / 1000
R = R_mm / 1000

# --- Step 2: Calculate intermediate terms for clarity ---
# Full permittivity of the insulator
epsilon = epsilon_r * epsilon_0
# Numerator of the main formula
numerator = 2 * math.pi * epsilon
# Term inside the logarithm
log_term_num = R**3 - m**3
log_term_den = 3 * r_w * m * R
log_term = log_term_num / log_term_den
# Denominator of the main formula
denominator = math.log(log_term)

# --- Step 3: Calculate capacitance in F/m ---
C_F_per_m = numerator / denominator

# --- Step 4: Convert to microfarads per km (uF/km) ---
# Conversion factor: 1 F/m = 1e6 uF / 1e-3 km = 1e9 uF/km
# Or step-by-step: (F * 1e6 uF/F) / (m * 1e-3 km/m)
C_uF_per_km = C_F_per_m * 1e6 * 1000

# --- Step 5: Print the final equation and result ---
print("The formula for capacitance per phase is: C = (2 * pi * epsilon) / ln((R^3 - m^3) / (3 * r_w * m * R))")
print("\n--- Calculation Steps ---")
print(f"1. Given values in SI units:")
print(f"   r_w = {r_w:.5f} m")
print(f"   m   = {m:.5f} m")
print(f"   R   = {R:.5f} m")
print(f"   epsilon = {epsilon_r} * {epsilon_0:.3e} = {epsilon:.5e} F/m")
print("\n2. The final equation with values plugged in:")
print(f"C = (2 * pi * {epsilon:.5e}) / ln(({R:.5f}^3 - {m:.5f}^3) / (3 * {r_w:.5f} * {m:.5f} * {R:.5f}))")
print("\n3. Evaluating the terms:")
print(f"C = ({numerator:.5e}) / ln(({log_term_num:.5e}) / ({log_term_den:.5e}))")
print(f"C = ({numerator:.5e}) / ln({log_term:.4f})")
print(f"C = ({numerator:.5e}) / {denominator:.4f}")
print(f"C = {C_F_per_m:.5e} F/m")
print("\n4. Converting to the final units (μF/km):")
print(f"C = {C_F_per_m:.5e} F/m * 1,000,000 μF/F * 1000 m/km")
print("\n--- Final Result ---")
print(f"The capacitance of the cable is {C_uF_per_km:.4f} μF/km.")
