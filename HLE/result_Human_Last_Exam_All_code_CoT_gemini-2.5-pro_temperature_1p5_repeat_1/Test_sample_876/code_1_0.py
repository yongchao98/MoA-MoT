import math

# --- Given values and constants ---
# Radius of the wires in meters
r_w = 11.25 / 1000
# Distance from wire center to cable center in meters
m = 17.32 / 1000
# External radius of the cable (inner radius of screen) in meters
R = 32.32 / 1000
# Relative permittivity of the insulator
epsilon_r = 4.2
# Permittivity of free space in F/m
epsilon_0 = 8.854e-12
# Pi
pi = math.pi

# --- Calculation ---
# 1. Calculate the numerator of the main formula
numerator = 2 * pi * epsilon_0 * epsilon_r

# 2. Calculate the components of the logarithm's argument
log_arg_numerator = R**3 - m**3
log_arg_denominator = 3 * r_w * m * R
log_argument = log_arg_numerator / log_arg_denominator

# 3. Calculate the denominator of the main formula (the natural log)
denominator = math.log(log_argument)

# 4. Calculate capacitance in Farads per meter (F/m)
C_F_per_m = numerator / denominator

# 5. Convert capacitance to microFarads per kilometer (uF/km)
C_uF_per_km = C_F_per_m * 1e9

# --- Output the results ---
print("This script calculates the capacitance of a three-phase cable with a common screen.")
print("\nThe formula used is: C = (2 * pi * epsilon_0 * epsilon_r) / ln((R^3 - m^3) / (3 * r_w * m * R))")
print("\n--- Final Equation with Substituted Values ---")
# The f-string below formats and prints the entire equation with all the numbers.
# We show the values in meters for physical consistency in the equation.
print(f"C = (2 * {pi:.5f} * {epsilon_0} * {epsilon_r}) / ln((({R})^3 - ({m})^3) / (3 * {r_w} * {m} * {R}))")
print(f"\nWhich simplifies to:")
print(f"C = ({numerator:.5e}) / ln({log_argument:.5f})")
print(f"C = ({numerator:.5e}) / {denominator:.5f}")
print(f"C = {C_F_per_m:.5e} F/m")
print("\n--- Final Answer ---")
print(f"The capacitance of the cable is: {C_uF_per_km:.4f} \u03BCF/km")
print("<<<0.5637>>>")