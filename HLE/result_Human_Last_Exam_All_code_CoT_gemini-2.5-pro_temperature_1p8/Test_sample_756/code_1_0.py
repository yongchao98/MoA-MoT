from fractions import Fraction

# The optimal location for the vertex of the parabola is x_0 = 1/2.
# This was found by parameterizing the coefficients a, b, c in terms of x_0
# and maximizing the objective function b+c.
x_0 = Fraction(1, 2)

# Based on the derivation, the coefficients b and c can be expressed in terms of x_0.
# b = 4*x_0 / (1 + x_0)**2
# c = (1 + 2*x_0 - x_0**2) / (1 + x_0)**2

# Calculate b and c for x_0 = 1/2
b = (4 * x_0) / (1 + x_0)**2
c = (1 + 2 * x_0 - x_0**2) / (1 + x_0)**2

# Calculate the value of |b| + |c|. Since we chose b>=0 and c>=0, this is b+c.
max_value = abs(b) + abs(c)

# The final equation is the sum of the absolute values of b and c.
# The problem asks to output each number in the final equation.
print(f"The analysis leads to optimal coefficients b = {b} and c = {c}.")
print(f"The maximum value of |b| + |c| is given by the sum:")
print(f"{abs(b)} + {abs(c)} = {max_value}")
