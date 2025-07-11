from fractions import Fraction

# Step 1: Define the coefficients of the candidate polynomial
# P(x) = ax^2 + bx + c = (-4/5)x^2 + (4/5)x + (4/5)
a = Fraction(-4, 5)
b = Fraction(4, 5)
c = Fraction(4, 5)

# Step 2: Verify that the polynomial satisfies the condition |P(x)| <= 1 for x in [-1, 1]
# This is a downward parabola, so we check the vertex and the endpoints.

# Vertex x-coordinate
x_vertex = -b / (2 * a)

# Value at the vertex (maximum value)
p_vertex = a * x_vertex**2 + b * x_vertex + c

# Values at the endpoints
p_at_1 = a * 1**2 + b * 1 + c
p_at_minus_1 = a * (-1)**2 + b * (-1) + c

# The minimum value on the interval is the minimum of the values at the endpoints
min_val = min(p_at_1, p_at_minus_1)
max_val = p_vertex

print("Analysis of the polynomial P(x) = (-4/5)x^2 + (4/5)x + (4/5) on [-1, 1]:")
print(f"Maximum value (at vertex x={x_vertex}): {max_val}")
print(f"Minimum value (at x={-1 if p_at_minus_1 < p_at_1 else 1}): {min_val}")

if max_val <= 1 and min_val >= -1:
    print("The condition |P(x)| <= 1 is satisfied.")
else:
    print("The condition |P(x)| <= 1 is NOT satisfied.")

# Step 3: Calculate the value of |b| + |c|
abs_b = abs(b)
abs_c = abs(c)
result = abs_b + abs_c

print("\nCalculating the maximum value of |b| + |c|:")
print(f"|b| + |c| = |{b}| + |{c}| = {abs_b} + {abs_c} = {result}")
