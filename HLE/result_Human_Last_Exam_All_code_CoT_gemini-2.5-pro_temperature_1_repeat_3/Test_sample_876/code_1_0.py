import math

# --- Given values and constants ---
r_w = 11.25      # Radius of the wires in mm
m = 17.32        # Distance from wire center to cable center in mm
R = 32.32        # External radius of the cable in mm
epsilon_r = 4.2  # Relative permittivity of the insulator
epsilon_0 = 8.854e-12 # Permittivity of free space in F/m

# --- Calculation Steps ---

# Calculate the argument of the natural logarithm
# The mm units cancel out, so we can use the values directly.
log_arg_num = R**3 - m**3
log_arg_den = 3 * r_w * m * R
log_arg = log_arg_num / log_arg_den

# Calculate the natural logarithm
ln_val = math.log(log_arg)

# Calculate the numerator of the capacitance formula, including the conversion factor
# The factor 1e9 converts from F/m to uF/km (1e6 for microfarads, 1e3 for km)
numerator_val = 2 * math.pi * epsilon_0 * epsilon_r * 1e9

# Calculate the final capacitance in uF/km
capacitance = numerator_val / ln_val

# --- Output the results ---
print("This script calculates the capacitance of a three-phase cable in uF/km.")
print("\nFinal Equation with values:")
# Print the full equation with all numbers substituted in
print(f"C = (2 * {math.pi:.5f} * {epsilon_0} * {epsilon_r} * 1e9) / ln(({R}**3 - {m}**3) / (3 * {r_w} * {m} * {R}))")
print(f"C = ({numerator_val:.2f}) / ln(({log_arg_num:.2f}) / ({log_arg_den:.2f}))")
print(f"C = {numerator_val:.2f} / ln({log_arg:.5f})")
print(f"C = {numerator_val:.2f} / {ln_val:.5f}")
print(f"C = {capacitance:.4f} uF/km")
<<<0.5655>>>