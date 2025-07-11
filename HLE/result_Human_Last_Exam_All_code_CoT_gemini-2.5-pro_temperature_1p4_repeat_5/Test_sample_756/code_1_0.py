from fractions import Fraction

# Based on the analytical derivation, one polynomial that maximizes |b|+|c|
# subject to the constraint is P(x) = -(8/9)x^2 + (8/9)x + 7/9.
# Let's define its coefficients.
a = -8.0/9.0
b = 8.0/9.0
c = 7.0/9.0

# The expression to maximize is |b| + |c|.
# Let's calculate the values of |b|, |c| and their sum.
abs_b = abs(b)
abs_c = abs(c)
max_value = abs_b + abs_c

# To verify, the range of P(x) on [-1, 1] must be within [-1, 1].
# For a downward parabola with vertex in (-1,1), the max is at the vertex
# and the min is at an endpoint.
x_vertex = -b / (2 * a)
p_max = a * x_vertex**2 + b * x_vertex + c
p_min = min(a - b + c, a + b + c) # P(-1) and P(1)

# We can assert that the found polynomial respects the constraints.
# Note: due to floating point arithmetic, we use a small tolerance.
assert -1.0 - 1e-9 <= p_min and p_max <= 1.0 + 1e-9

print("An optimal polynomial that maximizes |b|+|c| has been found.")
print(f"Its coefficients are a = {a:.4f}, b = {b:.4f}, c = {c:.4f}")
print("\nThe value of |b| + |c| for this polynomial gives the maximum value.")
print("The final equation is:")

# Presenting the result with fractions for exactness
b_frac = Fraction(b).limit_denominator()
c_frac = Fraction(c).limit_denominator()
max_val_frac = Fraction(max_value).limit_denominator()

print(f"|{b_frac}| + |{c_frac}| = {b_frac} + {c_frac} = {max_val_frac}")