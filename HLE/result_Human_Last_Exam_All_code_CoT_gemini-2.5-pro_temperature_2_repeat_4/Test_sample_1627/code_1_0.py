import sympy

# Define symbolic variables for inner radius 'a' and outer radius 'b'
a, b = sympy.symbols('a b')

# Define the relationship between the outer and inner radius
# b = 2a
radius_relationship = {b: 2*a}

# The formula for the ratio of maximum tangential stress to internal pressure is:
# Ratio = (b^2 + a^2) / (b^2 - a^2)
numerator_expr = b**2 + a**2
denominator_expr = b**2 - a**2

# Substitute the radius relationship into the expressions
final_numerator = numerator_expr.subs(radius_relationship)
final_denominator = denominator_expr.subs(radius_relationship)

# Simplify the expressions to remove the variable 'a'
# final_numerator becomes 5*a**2, so the numeric part is 5
# final_denominator becomes 3*a**2, so the numeric part is 3
num_val = sympy.simplify(final_numerator / a**2)
den_val = sympy.simplify(final_denominator / a**2)

# Calculate the final decimal value of the ratio
ratio_value = num_val / den_val

print("The formula for the ratio of maximum tangential stress to internal pressure is (b^2 + a^2) / (b^2 - a^2).")
print("Given that the outer radius 'b' is twice the inner radius 'a' (b = 2a), we substitute this into the formula.")
print("Ratio = ((2a)^2 + a^2) / ((2a)^2 - a^2) = (5a^2) / (3a^2)")
print(f"After simplifying, the final ratio equation is:")
# The final prompt asked to output each number in the final equation.
# So we print the numerator and denominator separately.
print(f"Maximum Tangential Stress / Internal Pressure = {num_val} / {den_val}")
print(f"The decimal value of this ratio is approximately {float(ratio_value):.3f}")
