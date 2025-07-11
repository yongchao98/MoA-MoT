import math

# Description:
# This script calculates the capacitance of a three-phase cable with a common screen.
# The result is provided in microfarads per kilometer (uF/km).

# --- Given values and constants ---
r_w = 11.25         # Wire radius in mm
m = 17.32           # Distance from wire center to cable center in mm
R = 32.32           # External radius of the cable in mm
epsilon_r = 4.2     # Relative permittivity of the insulator (dimensionless)
epsilon_0 = 8.854e-12 # Permittivity of free space in F/m
pi = math.pi

# --- Calculation ---
# The units (mm) cancel out in the logarithm argument.
num_log_arg = R**3 - m**3
den_log_arg = 3 * m * r_w * (R - m)
log_arg = num_log_arg / den_log_arg
ln_val = math.log(log_arg)

# Numerator of the main capacitance formula
num_C_formula = 2 * pi * epsilon_0 * epsilon_r

# Capacitance in F/m
C_F_per_m = num_C_formula / ln_val

# Convert from F/m to uF/km by multiplying by 10^9
conversion_factor = 1e9
C_uF_per_km = C_F_per_m * conversion_factor

# --- Output ---
print("The equation for capacitance (C) in uF/km is:")
print("C = (2 * pi * epsilon_0 * epsilon_r) / ln((R^3 - m^3) / (3 * m * r_w * (R - m))) * 1e9\n")

print("Substituting the given values into the equation:")
print(f"C = (2 * {pi:.6f} * {epsilon_0} * {epsilon_r}) / ln(({R}**3 - {m}**3) / (3 * {m} * {r_w} * ({R} - {m}))) * {conversion_factor:.0e}")
print("   Solving the equation yields:")
print(f"C = ({num_C_formula:.6e}) / ln({log_arg:.6f}) * {conversion_factor:.0e}")
print(f"C = ({num_C_formula:.6e}) / {ln_val:.6f} * {conversion_factor:.0e}")
print(f"C = {C_uF_per_km:.4f} \u03BCF/km\n")
print(f"The capacitance of the cable is {C_uF_per_km:.4f} \u03BCF/km.")
<<<0.1981>>>