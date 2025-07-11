import sympy

# Define the symbols for the square roots
sqrt2 = sympy.sqrt(2)
sqrt3 = sympy.sqrt(3)
sqrt6 = sympy.sqrt(6)

# The expressions from the commutativity check
# LHS: (3-sqrt(3)) / ((2+sqrt(2))*sqrt(3))
lhs = (3 - sqrt3) / ((2 + sqrt2) * sqrt3)

# RHS: (2-sqrt(2))*sqrt(3) / (3+sqrt(3))
rhs = (2 - sqrt2) * sqrt3 / (3 + sqrt3)

# Simplify both expressions to a canonical form
simplified_lhs = sympy.simplify(lhs)
simplified_rhs = sympy.simplify(rhs)

# The equation from cross-multiplication
# LHS: (3-sqrt(3))*(3+sqrt(3))
eq_lhs = (3 - sqrt3) * (3 + sqrt3)
# RHS: (2-sqrt(2))*sqrt(3)*(2+sqrt(2))*sqrt(3)
eq_rhs = (2 - sqrt2) * sqrt3 * (2 + sqrt2) * sqrt3

print(f"Let a = (2-sqrt(2))*sqrt(3)*(2+sqrt(2)) and b = (sqrt(3)-1)*(3+sqrt(3))")
# I use the identity derived from cross multiplying without denominator: (2-sqrt(2))sqrt(3)(2+sqrt(2)) = (sqrt(3)-1)(3+sqrt(3))
# which is equivalent to 2*sqrt(3) = 2*sqrt(3)
# In my thoughts: ((2-sqrt(2))*sqrt(3)*(2+sqrt(2))) - ((sqrt(3)-1)*(3+sqrt(3))) must be zero.
# Left side of equation
left_val = (2-sqrt2)*sqrt3*(2+sqrt2)
# Right side of equation
right_val = (sqrt3-1)*(3+sqrt3)
print(f"The equation we test is: {sympy.simplify(left_val)} = {sympy.simplify(right_val)}")
# The final numbers are simplified. Let's show the equation holds by simplifying their difference to zero.
difference = sympy.simplify(left_val - right_val)
print(f"The simplified difference is: {difference}")
if difference == 0:
    print("The identity holds, so the generators commute.")
else:
    print("The identity does not hold.")
