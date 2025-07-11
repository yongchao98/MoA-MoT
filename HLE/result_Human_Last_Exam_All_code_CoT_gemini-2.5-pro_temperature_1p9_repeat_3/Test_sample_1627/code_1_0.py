import sympy

# Define symbolic variables for clarity
a, b, Pi = sympy.symbols('a b P_i')

# The relationship between outer and inner radius
# b = 2*a
b_val = 2
a_val = 1

# Numerator and denominator of the ratio formula
numerator_expr = b**2 + a**2
denominator_expr = b**2 - a**2

# Substitute the values
numerator_val = numerator_expr.subs({b: b_val, a: a_val})
denominator_val = denominator_expr.subs({b: b_val, a: a_val})

# The ratio is the division of the two
ratio_val = sympy.Rational(numerator_val, denominator_val)

# Print the final equation and the result
print("The formula for the ratio of maximum tangential stress to internal pressure is:")
print("Ratio = (b^2 + a^2) / (b^2 - a^2)")
print("\nGiven the outer radius (b) is twice the inner radius (a), let's set a = 1 and b = 2.")
print("\nSubstituting the values into the formula:")
# Note: f-strings in Python don't directly handle the superscript '²'. Using '**2' for clarity in the code representation.
final_equation = f"Ratio = ({b_val}**2 + {a_val}**2) / ({b_val}**2 - {a_val}**2)"
final_calculation = f"Ratio = ({b_val**2} + {a_val**2}) / ({b_val**2} - {a_val**2})"
final_result = f"Ratio = {numerator_val} / {denominator_val}"
print(final_equation)
print(final_calculation)
print(final_result)
print(f"\nThus, the ratio is {ratio_val}, which is approximately {float(ratio_val):.3f}.")