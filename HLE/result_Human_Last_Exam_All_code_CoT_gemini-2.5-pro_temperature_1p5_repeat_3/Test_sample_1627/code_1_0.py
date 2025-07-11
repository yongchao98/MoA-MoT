import sympy

# Define symbolic variables
a = sympy.Symbol('a')
b = sympy.Symbol('b')

# Define the ratio formula based on Lame's equation
# Ratio = (sigma_t_max / P_i) = (b^2 + a^2) / (b^2 - a^2)
numerator_expr = b**2 + a**2
denominator_expr = b**2 - a**2
ratio_expr = numerator_expr / denominator_expr

# Apply the given condition: outer radius is twice the inner radius (b = 2a)
ratio_with_condition = ratio_expr.subs(b, 2*a)

# Simplify the expression to find the numerical ratio
final_ratio = sympy.simplify(ratio_with_condition)

# To display the calculation step-by-step, let's substitute numerical values
# We can set a=1 for simplicity, which makes b=2
a_val = 1
b_val = 2

# Calculate the numerator and denominator with these values
num_val = b_val**2 + a_val**2
den_val = b_val**2 - a_val**2

print("In the context of elasticity theory, the ratio of maximum tangential stress to internal pressure is given by:")
print("Ratio = (b^2 + a^2) / (b^2 - a^2)")
print("where 'a' is the inner radius and 'b' is the outer radius.")
print("\nGiven that the outer radius is twice the inner radius, we have b = 2a.")
print("For calculation purposes, we can set a = 1 and b = 2.")
print("\nSubstituting these values into the formula:")
# Note: The final print statement is formatted to show each number explicitly as requested.
print(f"Ratio = ({b_val}^2 + {a_val}^2) / ({b_val}^2 - {a_val}^2)")
print(f"Ratio = ({b_val**2} + {a_val**2}) / ({b_val**2} - {a_val**2})")
print(f"Ratio = {num_val} / {den_val}")

# The final_ratio variable from sympy gives the exact fraction
print(f"\nThe simplified ratio is {final_ratio}.")
# Also print the decimal representation
print(f"As a decimal, the ratio is approximately {float(final_ratio):.3f}.")