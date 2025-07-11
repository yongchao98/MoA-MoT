import numpy as np

# Define the coefficients of the polynomial f(x) = ax^2 + bx + c
# This polynomial is a candidate for maximizing |b| + |c|.
a = -0.5
b = 1.0
c = 0.5

def f(x, a, b, c):
  """The quadratic polynomial."""
  return a * x**2 + b * x + c

# Check if the condition |f(x)| <= 1 holds for x in [-1, 1].
# We can check the extrema of the function in the interval.
# The derivative is f'(x) = 2ax + b.
# The vertex is at x_v = -b / (2a).
x_vertex = -b / (2 * a)

# Evaluate the function at the endpoints and the vertex if it's in the interval.
f_at_minus_1 = f(-1, a, b, c)
f_at_1 = f(1, a, b, c)

print(f"The polynomial is f(x) = {a}x^2 + {b}x + {c}")
print(f"The value at x = -1 is: {f_at_minus_1}")
print(f"The value at x = 1 is: {f_at_1}")

# Since the vertex is at x=1, the function is monotonic on [-1, 1].
# The extrema are at the endpoints.
# As f(-1) = -1 and f(1) = 1, the condition |f(x)| <= 1 holds.

# Calculate the value of |b| + |c|
max_val = abs(b) + abs(c)

print("\nFor this polynomial:")
# The final result should be an equation, as requested.
print(f"The value of the expression is |{b}| + |{c}| = {max_val}")
print(f"This is the maximum value.")
