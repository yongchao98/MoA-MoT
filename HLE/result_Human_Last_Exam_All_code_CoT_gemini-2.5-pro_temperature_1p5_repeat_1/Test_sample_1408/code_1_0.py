import math

# The maximal overhang for three cubes of side length 1 is 11*d/6,
# where d is the distance from the center to a corner of a cube's face.
# For a cube of side length 1, d = 1/sqrt(2).
# Overhang = 11 / (6 * sqrt(2)) = (11 * sqrt(2)) / 12.

# We need to express this in the form (a + sqrt(b)) / (1 + c).

# Denominator is 12, so 1 + c = 12.
denominator = 12
c = denominator - 1

# Numerator is 11 * sqrt(2). We write this as sqrt(b).
# a + sqrt(b) = 11 * sqrt(2)
# Since a must be a non-negative integer, a must be 0.
a = 0
# So, sqrt(b) = 11 * sqrt(2) = sqrt(11^2 * 2) = sqrt(242).
b = 11**2 * 2

print("The determined non-negative integers are:")
print(f"a = {a}")
print(f"b = {b}")
print(f"c = {c}")

print("\nThe final equation for the maximal overhang is:")
# Using 'math.sqrt' for the numerical value, and integer 'b' for the expression.
print(f"({a} + sqrt({b})) / (1 + {c}) = {(a + math.sqrt(b)) / (1 + c):.4f}")
