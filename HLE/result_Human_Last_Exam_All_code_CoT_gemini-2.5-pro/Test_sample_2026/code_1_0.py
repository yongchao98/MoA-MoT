import math

# The problem is reduced to finding the surface area of a sphere defined by the values of (x0, y0, z0).
# The derivation shows that these values satisfy the equation:
# x0^2 + y0^2 + z0^2 = R^2
# where R^2 = 10^25 / 2.
# The surface area of this sphere is given by the formula A = 4 * pi * R^2.
# Substituting R^2, we get A = 4 * pi * (10^25 / 2) = 2 * pi * 10^25.

# Constants from the simplified final equation for the area: Area = factor * pi * base^exponent
factor = 2
base = 10
exponent = 25
power_of_ten = base**exponent

# Calculate the final area
area = factor * math.pi * power_of_ten

# As requested, we print each number in the final equation.
# The final equation is: Area = 2 * pi * 10^25
print("The final equation for the area is of the form: A = (factor) * pi * (power_of_ten)")
print(f"The factor is: {factor}")
print(f"Pi is approximately: {math.pi}")
print(f"The power of ten is: {power_of_ten:.0e}")
print("\n---")
print(f"The calculated area is: {area:.12e}")