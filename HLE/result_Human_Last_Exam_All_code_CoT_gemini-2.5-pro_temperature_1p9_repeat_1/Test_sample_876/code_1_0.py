import math

# --- Task: Calculate Capacitance of a Three-Phase Cable ---

# 1. Define physical constants and given parameters
epsilon_0 = 8.854e-12  # Permittivity of free space in F/m

# Parameters from the problem statement
r_w_mm = 11.25    # Wire radius in mm
m_mm = 17.32      # Distance from wire center to cable center in mm
R_mm = 32.32      # External radius of the cable in mm
epsilon_r = 4.2   # Relative permittivity of the insulator

# 2. Convert all dimensions to SI units (meters)
r_w = r_w_mm / 1000.0
m = m_mm / 1000.0
R = R_mm / 1000.0

# 3. Calculate the capacitance in Farads per meter (F/m)
# The formula for per-phase capacitance is:
# C = (2 * pi * epsilon_r * epsilon_0) / ln( (R^3 - m^3) / (3 * R * m * r_w) )
epsilon = epsilon_r * epsilon_0
numerator = 2 * math.pi * epsilon
denominator_arg = (R**3 - m**3) / (3 * R * m * r_w)
denominator = math.log(denominator_arg)
C_F_per_m = numerator / denominator

# 4. Convert the result to microfarads per kilometer (uF/km)
# Conversion factor: 1 F/m = 10^9 uF/km
C_uF_per_km = C_F_per_m * 1e9

# 5. Output the results, including the equation with values
print("--- Calculation Details ---")
print("The formula for capacitance per phase (C) is:")
print("C = (2 * pi * epsilon_r * epsilon_0) / ln((R^3 - m^3) / (3 * R * m * r_w))\n")

print("The final equation with the given values (in SI units) is:")
# The following print statement fulfills the requirement to show each number in the equation.
print(f"C (F/m) = (2 * pi * {epsilon_r} * {epsilon_0:.4e}) / ln((({R:.5f})^3 - ({m:.5f})^3) / (3 * {R:.5f} * {m:.5f} * {r_w:.5f}))\n")

print("--- Result ---")
print(f"The capacitance of the cable is {C_uF_per_km:.4f} \u03BCF/km.")
<<<0.5646>>>