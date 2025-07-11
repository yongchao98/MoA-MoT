import numpy as np

# The problem asks for the determinant of the matrix A.
# The information "Markov quiver for m_{1/2} = 13" is context for where this matrix arises,
# but it is not needed for the computation itself.

# The given matrix is:
A = np.array([
    [2, -1, -1],
    [-1, 2, -1],
    [-1, -1, 2]
])

# For a 3x3 matrix [[a, b, c], [d, e, f], [g, h, i]], the determinant is calculated as:
# det = a(ei - fh) - b(di - fg) + c(dh - eg)

# Let's extract the elements from matrix A to use in the formula.
a, b, c = A[0, 0], A[0, 1], A[0, 2]
d, e, f = A[1, 0], A[1, 1], A[1, 2]
g, h, i = A[2, 0], A[2, 1], A[2, 2]

print(f"The given matrix A is:")
print(A)
print("\nWe calculate the determinant using the formula: a(ei - fh) - b(di - fg) + c(dh - eg)")

# Step 1: Substitute the numbers from matrix A into the formula.
# Note that the formula has a '-b' term, which becomes '-(-1)' = '+1'.
print("\nStep 1: Substituting the numbers into the formula gives:")
print(f"Determinant = ({a}) * (({e})*({i}) - ({f})*({h})) - ({b}) * (({d})*({i}) - ({f})*({g})) + ({c}) * (({d})*({h}) - ({e})*({g}))")

# Step 2: Calculate the values inside each parenthesis (the cofactors).
val1 = e * i - f * h
val2 = d * i - f * g
val3 = d * h - e * g

print("\nStep 2: Evaluating the expressions inside the parentheses:")
print(f"({e})*({i}) - ({f})*({h}) = {4 - 1} = {val1}")
print(f"({d})*({i}) - ({f})*({g}) = {-2 - 1} = {val2}")
print(f"({d})*({h}) - ({e})*({g}) = {1 - (-2)} = {val3}")

# Step 3: Substitute these values back into the equation.
print("\nStep 3: The equation becomes:")
# To correctly show the sign for the 'b' term, we calculate -b
term_b_sign = -b
# Display it nicely with a plus sign if positive
sign = '+' if term_b_sign > 0 else ''
print(f"Determinant = ({a})*({val1}) {sign} ({term_b_sign})*({val2}) + ({c})*({val3})")

# Step 4: Perform the final multiplications and sum.
final_result = a * val1 - b * val2 + c * val3
print("\nStep 4: Performing the final calculation:")
print(f"Determinant = {a * val1} + {term_b_sign * val2} + {c * val3}")
print(f"Determinant = {6} - {3} - {3} = {final_result}")

print("\n-------------------------------------------")
print(f"The determinant of the adjacency matrix is: {final_result}")
print("-------------------------------------------")