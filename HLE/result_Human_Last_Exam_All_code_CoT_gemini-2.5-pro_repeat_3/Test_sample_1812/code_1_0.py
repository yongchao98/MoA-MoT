import numpy as np

# The matrix provided in the problem
A = np.array([
    [2, -1, -1],
    [-1, 2, -1],
    [-1, -1, 2]
])

# Extracting elements of the matrix for the formula
# a b c
# d e f
# g h i
a, b, c = A[0, 0], A[0, 1], A[0, 2]
d, e, f = A[1, 0], A[1, 1], A[1, 2]
g, h, i = A[2, 0], A[2, 1], A[2, 2]

# Calculate the determinant using the standard formula for a 3x3 matrix
determinant_value = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)

# Print the calculation steps with all numbers shown in the equation
print("To compute the determinant of the matrix A:")
print(A)
print("\nWe use the formula for a 3x3 matrix: det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)")
print("\nSubstituting the values from matrix A into the formula:")
print(f"det(A) = ({a}) * (({e})*({i}) - ({f})*({h})) - ({b}) * (({d})*({i}) - ({f})*({g})) + ({c}) * (({d})*({h}) - ({e})*({g}))")
print(f"       = ({a}) * ({e*i} - {f*h}) - ({b}) * ({d*i} - {f*g}) + ({c}) * ({d*h} - {e*g})")
print(f"       = ({a}) * ({e*i - f*h}) - ({b}) * ({d*i - f*g}) + ({c}) * ({d*h - e*g})")
print(f"       = {a * (e*i - f*h)} - {b * (d*i - f*g)} + {c * (d*h - e*g)}")
print(f"       = {determinant_value}")

print(f"\nThe final determinant of the matrix is: {determinant_value}")