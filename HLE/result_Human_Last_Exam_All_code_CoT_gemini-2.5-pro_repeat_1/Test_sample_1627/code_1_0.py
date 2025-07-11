import sympy

# Define symbols for the radii and pressure
a, b, P_i = sympy.symbols('a b P_i')

# 1. Formula for maximum tangential stress (sigma_t_max)
# sigma_t_max = P_i * (b**2 + a**2) / (b**2 - a**2)

# 2. The required ratio is sigma_t_max / P_i
ratio_expr = (b**2 + a**2) / (b**2 - a**2)

# 3. Apply the condition that the outer radius is twice the inner radius (b = 2a)
# For calculation, we can set a=1 and b=2, as the 'a' term will cancel out.
inner_radius_val = 1
outer_radius_val = 2

# Substitute the values into the expression
# Numerator calculation
num_val = outer_radius_val**2 + inner_radius_val**2
# Denominator calculation
den_val = outer_radius_val**2 - inner_radius_val**2

# Calculate the final ratio
final_ratio = num_val / den_val

# 4. Print the steps of the calculation
print("The ratio of maximum tangential stress (σ_t_max) to internal pressure (P_i) is given by the formula:")
print("Ratio = (b² + a²) / (b² - a²)")
print("\nGiven that the outer radius (b) is twice the inner radius (a), we can set a=1 and b=2.")
print("\nSubstituting these values into the formula:")
print(f"Ratio = ({outer_radius_val}² + {inner_radius_val}²) / ({outer_radius_val}² - {inner_radius_val}²)")
print(f"Ratio = ({outer_radius_val**2} + {inner_radius_val**2}) / ({outer_radius_val**2} - {inner_radius_val**2})")
print(f"Ratio = {num_val} / {den_val}")
print(f"\nTherefore, the final ratio is {final_ratio:.3f}, or as a fraction, 5/3.")
