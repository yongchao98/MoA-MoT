import math

# --- Given values ---
r_w_mm = 11.25  # Wire radius in mm
m_mm = 17.32    # Distance from wire center to cable center in mm
R_mm = 32.32    # External radius of cable insulation in mm
epsilon_r = 4.2 # Relative permittivity of the insulator

# --- Constants ---
epsilon_0 = 8.854e-12 # Permittivity of free space in F/m

# --- Calculations ---
# Step 1: Convert all lengths from mm to m for unit consistency
r_w = r_w_mm / 1000.0
m = m_mm / 1000.0
R = R_mm / 1000.0

# Step 2: Calculate the components of the formula
# Numerator of the main formula
numerator = 2 * math.pi * epsilon_0 * epsilon_r

# Argument of the natural logarithm
ln_arg_num = R**2 - m**2
ln_arg_den = r_w * R
ln_arg = ln_arg_num / ln_arg_den

# Denominator of the main formula (natural logarithm)
denominator = math.log(ln_arg)

# Capacitance in Farads per meter (F/m)
C_F_per_m = numerator / denominator

# Step 3: Convert the final result to microfarads per kilometer (uF/km)
# 1 F/m = 10^9 uF/km
C_uF_per_km = C_F_per_m * 1e9

# --- Output ---
print("Calculation for the capacitance of a three-phase cable.")
print("The formula used is: C = (2 * pi * epsilon_0 * epsilon_r) / ln((R^2 - m^2) / (r_w * R))")
print("\n--- The final equation with the given values ---")
# Print each number used in the final equation
print("C (uF/km) = (2 * {} * {:.3e} F/m * {}) / ln(({:.5f} m)^2 - ({:.5f} m)^2) / (({:.5f} m) * ({:.5f} m)) * 1e9".format(
    round(math.pi, 4), epsilon_0, epsilon_r, R, m, r_w, R
))
print("\n--- Final Result ---")
print(f"The capacitance of the cable is {C_uF_per_km:.3f} uF/km.")
<<<0.326>>>