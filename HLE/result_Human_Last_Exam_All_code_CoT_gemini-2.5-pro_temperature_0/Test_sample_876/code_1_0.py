import math

# --- Given values ---
r_w = 11.25  # wire radius in mm
m = 17.32    # distance from center of a wire to the center of the cable in mm
R = 32.32    # external radius of the cable in mm
epsilon_r = 4.2 # relative permittivity of the insulator

# --- Constants ---
epsilon_0 = 8.854e-12 # permittivity of free space in F/m
pi = math.pi

# --- Calculation ---
# The geometric part of the formula is a ratio, so we can use mm directly.
log_numerator = R**3 - m**3
log_denominator = 3 * m * R * r_w
log_argument = log_numerator / log_denominator

# Denominator of the main formula
final_denominator = math.log(log_argument)

# Numerator of the main formula, including the conversion factor to uF/km
# Conversion factor from F/m to uF/km is 10^9
final_numerator = 2 * pi * epsilon_r * epsilon_0 * 1e9

# Final capacitance in uF/km
C_uF_per_km = final_numerator / final_denominator

# --- Output ---
print("To find the capacitance C in microfarads per kilometer (uF/km), we use the formula:")
print("C = [ (2 * pi * epsilon_r * epsilon_0) / ln((R^3 - m^3) / (3 * m * R * r_w)) ] * 10^9\n")

print("Substituting the given values (with lengths in mm):")
print(f"C = [ (2 * {pi:.5f} * {epsilon_r} * {epsilon_0}) / ln(({R}^3 - {m}^3) / (3 * {m} * {R} * {r_w})) ] * 10^9")

print("\nCalculating the numerator and denominator:")
print(f"C = {final_numerator:.5f} / {final_denominator:.5f}")

print("\nThe final capacitance is:")
print(f"C = {C_uF_per_km:.4f} uF/km")