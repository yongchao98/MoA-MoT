import math

# Given values
m1 = 1
m2 = 2
R = 3
d = 1
g = 10

# The task is to propose an equation for v0. We will construct the string for the equation
# using the derived formula:
# v0 = sqrt[ (4 * g * (m1*R + m2*d) * (m1*R^2 + m2*d^2)) / (m2^2 * d^2) ]
# And we will print the equation with the numbers substituted in.

# The prompt asks for an equation for v0, not v0^2, so we'll use the form with the square root.
# Let's write the equation clearly showing each number.

equation_str = f"v0 = sqrt[(4 * {g} * ({m1} * {R} + {m2} * {d}) * ({m1} * {R}^2 + {m2} * {d}^2)) / ({m2}^2 * {d}^2)]"

print("The equation for v0 is:")
print(equation_str)

# Let's also calculate the numerical value for completeness, although it's not explicitly asked for.
term1 = m1*R + m2*d
term2 = m1*R**2 + m2*d**2
numerator = 4 * g * term1 * term2
denominator = m2**2 * d**2
v0_squared = numerator / denominator
v0 = math.sqrt(v0_squared)

# print("\nCalculation:")
# print(f"v0 = sqrt[(4 * {g} * {term1} * {term2}) / {denominator}]")
# print(f"v0 = sqrt[{numerator / denominator}]")
# print(f"v0 = {v0}")