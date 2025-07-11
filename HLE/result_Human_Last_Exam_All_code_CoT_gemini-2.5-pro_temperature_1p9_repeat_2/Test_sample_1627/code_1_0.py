import sympy

# Define the symbols
# P_i: Internal Pressure
# a: Inner radius
# b: Outer radius
# ratio: The ratio of maximum tangential stress to internal pressure
a, b, P_i = sympy.symbols('a b P_i')
sigma_t_max_over_Pi = (b**2 + a**2) / (b**2 - a**2)

# Apply the given condition: b = 2a
# We create a new expression by substituting b with 2*a
ratio_expr = sigma_t_max_over_Pi.subs(b, 2*a)

# Simplify the expression to find the numerical ratio
final_ratio = sympy.simplify(ratio_expr)

# Let's print the step-by-step derivation
print("The formula for the ratio of maximum tangential stress to internal pressure is:")
print("(σ_t_max / P_i) = (b^2 + a^2) / (b^2 - a^2)")
print("\nGiven the condition that the outer radius is twice the inner radius (b = 2a), we substitute this into the formula:")
# Get the numerator and denominator from the substituted expression before simplification
numerator_expr = sympy.denom(ratio_expr)
denominator_expr = sympy.numer(ratio_expr)
print(f"(σ_t_max / P_i) = ( (2a)^2 + a^2 ) / ( (2a)^2 - a^2 )")
print(f"(σ_t_max / P_i) = ( {sympy.expand(denominator_expr)} ) / ( {sympy.expand(numerator_expr)} )")

simplified_numerator = sympy.simplify(denominator_expr)
simplified_denominator = sympy.simplify(numerator_expr)
print(f"(σ_t_max / P_i) = ( {simplified_numerator} ) / ( {simplified_denominator} )")

# Cancel out the a^2 term
print(f"(σ_t_max / P_i) = {int(simplified_numerator/a**2)} / {int(simplified_denominator/a**2)}")

# Display the final ratio as a fraction and a decimal
print(f"\nThe final ratio is {final_ratio}.")
print(f"As a decimal, this is approximately {float(final_ratio):.3f}.")
print("\nThe final equation is:")
print(f"(4*a^2 + a^2) / (4*a^2 - a^2) = (5*a^2) / (3*a^2) = {int(simplified_numerator/a**2)} / {int(simplified_denominator/a**2)}")
