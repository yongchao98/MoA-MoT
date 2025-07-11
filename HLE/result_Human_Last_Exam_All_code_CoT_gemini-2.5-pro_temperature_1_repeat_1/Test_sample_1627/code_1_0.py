import sympy

# Define the symbols for the radii
a, b = sympy.symbols('a b')

# Define the relationship between the radii
# The outer radius (b) is twice the inner radius (a)
radius_relation = {b: 2*a}

# The formula for the ratio of maximum tangential stress to internal pressure
ratio_formula = (b**2 + a**2) / (b**2 - a**2)

# Substitute the radius relationship into the formula
# The 'a' terms will cancel out, giving a numerical result
final_ratio = ratio_formula.subs(radius_relation)

# Extract the numerator and denominator to display the calculation
numerator_expr = (b**2 + a**2).subs(radius_relation)
denominator_expr = (b**2 - a**2).subs(radius_relation)

# To show the numbers, we can substitute a=1
numerator_val = numerator_expr.subs({a: 1})
denominator_val = denominator_expr.subs({a: 1})

# Print the calculation steps
print("The formula for the ratio of maximum tangential stress to internal pressure is: (b^2 + a^2) / (b^2 - a^2)")
print("Given that the outer radius 'b' is twice the inner radius 'a' (b = 2a).")
print("Substituting b = 2a into the formula gives: ((2a)^2 + a^2) / ((2a)^2 - a^2)")
print(f"This simplifies to: ({numerator_expr}) / ({denominator_expr})")
print(f"The 'a^2' terms cancel, leaving the ratio. To show the numerical calculation, let's set a=1 and b=2:")
print(f"Ratio = (2^2 + 1^2) / (2^2 - 1^2)")
print(f"Ratio = (4 + 1) / (4 - 1)")
print(f"Ratio = {int(numerator_val)} / {int(denominator_val)}")

# The final answer is the fraction 5/3
print(f"\nThe final ratio is 5/3, which is approximately {float(final_ratio):.3f}.")
