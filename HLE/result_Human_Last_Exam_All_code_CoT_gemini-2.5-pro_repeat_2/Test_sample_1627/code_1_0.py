import sympy

# Let's define the symbols for our calculation
a, b, P_i = sympy.symbols('a b P_i')

# The formula for the ratio of maximum tangential stress to internal pressure is:
# Ratio = (b^2 + a^2) / (b^2 - a^2)
# We are given the condition that the outer radius is twice the inner radius.
condition = {b: 2*a}

# The numerator of the ratio formula
numerator_expr = b**2 + a**2
# The denominator of the ratio formula
denominator_expr = b**2 - a**2

# Substitute the condition b = 2a into the expressions
final_numerator = sympy.simplify(numerator_expr.subs(condition))
final_denominator = sympy.simplify(denominator_expr.subs(condition))

# Extract the numerical coefficients
num_coeff = final_numerator.as_coeff_Mul()[0]
den_coeff = final_denominator.as_coeff_Mul()[0]

print("The formula for the ratio of maximum tangential stress to internal pressure is:")
print("Ratio = (b^2 + a^2) / (b^2 - a^2)\n")

print("Given b = 2a, we substitute this into the formula:")
print("Ratio = ((2a)^2 + a^2) / ((2a)^2 - a^2)")
print("Ratio = (4a^2 + a^2) / (4a^2 - a^2)")
print("Ratio = 5a^2 / 3a^2\n")

print("After canceling the 'a^2' term, we get the final ratio.")
print("The final equation for the ratio is:")
print(f"Maximum Tangential Stress / Internal Pressure = {int(num_coeff)} / {int(den_coeff)}")
