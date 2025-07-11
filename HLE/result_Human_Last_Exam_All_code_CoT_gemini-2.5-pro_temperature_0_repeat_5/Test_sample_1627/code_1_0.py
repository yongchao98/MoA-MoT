import sympy

# Define the symbols for the radii
a, b = sympy.symbols('a b')

# Define the relationship between the outer and inner radius
# b = 2a
radius_relation = sympy.Eq(b, 2*a)

# The formula for the ratio of maximum tangential stress to internal pressure
ratio_formula = (b**2 + a**2) / (b**2 - a**2)

# Substitute the radius relationship into the formula
final_ratio = ratio_formula.subs(b, 2*a)

# Simplify the expression to get the numerical value
simplified_ratio = sympy.simplify(final_ratio)

# Print the calculation steps
print("The formula for the ratio of maximum tangential stress to internal pressure is: (b^2 + a^2) / (b^2 - a^2)")
print("Given that the outer radius 'b' is twice the inner radius 'a', we substitute b = 2a.")
print(f"Ratio = ((2a)^2 + a^2) / ((2a)^2 - a^2)")
print(f"      = (4a^2 + a^2) / (4a^2 - a^2)")
print(f"      = 5a^2 / 3a^2")
print(f"      = {simplified_ratio}")

# To show the decimal value
print(f"\nThe numerical value of the ratio is {float(simplified_ratio):.4f}")
