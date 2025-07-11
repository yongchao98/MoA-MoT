import sympy

# Define symbolic variables for the radii
a, b = sympy.symbols('a b')

# Define the relationship between the radii
# b = 2a
radius_relation = {b: 2*a}

# The formula for the ratio of maximum tangential stress to internal pressure is:
# (a^2 + b^2) / (b^2 - a^2)
numerator_expr = a**2 + b**2
denominator_expr = b**2 - a**2
ratio_expr = numerator_expr / denominator_expr

# Substitute the radius relation into the formula
ratio_substituted = ratio_expr.subs(radius_relation)

# Simplify the expression
simplified_ratio = sympy.simplify(ratio_substituted)

# To show the step-by-step calculation with numbers, let's substitute b=2 and a=1
b_val = 2
a_val = 1
numerator_val = b_val**2 + a_val**2
denominator_val = b_val**2 - a_val**2

print("The formula for the ratio of maximum tangential stress to internal pressure is:")
print("(σ_t_max) / P_i = (a² + b²) / (b² - a²)")
print("\nGiven the outer radius 'b' is twice the inner radius 'a', we set b = 2a.")
print("Substituting b = 2a into the formula gives:")
print("(σ_t_max) / P_i = (a² + (2a)²) / ((2a)² - a²)")
print("(σ_t_max) / P_i = (a² + 4a²) / (4a² - a²)")
print("(σ_t_max) / P_i = 5a² / 3a²")
print(f"(σ_t_max) / P_i = {simplified_ratio}")

print("\nFinal calculation:")
# Print the final equation with the numbers plugged in
print(f"The equation evaluates to: ({b_val**2} + {a_val**2}) / ({b_val**2} - {a_val**2})")
print(f"Which simplifies to: {numerator_val} / {denominator_val}")
print(f"The ratio is {float(simplified_ratio):.4f}")