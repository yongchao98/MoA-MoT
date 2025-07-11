import numpy as np

# The problem asks to compute the determinant of the adjacency matrix A.
# The matrix is given as:
A = np.array([
    [2, -1, -1],
    [-1, 2, -1],
    [-1, -1, 2]
])

# We can compute the determinant using the formula for a 3x3 matrix:
# det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)
# where the matrix is:
# [[a, b, c],
#  [d, e, f],
#  [g, h, i]]

a, b, c = A[0, 0], A[0, 1], A[0, 2]
d, e, f = A[1, 0], A[1, 1], A[1, 2]
g, h, i = A[2, 0], A[2, 1], A[2, 2]

# Calculate the determinant
determinant = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)

# Print the final equation with all the numbers
print("The determinant calculation is as follows:")
print(f"{a} * (({e} * {i}) - ({f} * {h})) - ({b}) * (({d} * {i}) - ({f} * {g})) + ({c}) * (({d} * {h}) - ({e} * {g}))")
print(f"= {a} * ({e * i} - {f * h}) - ({b}) * ({d * i} - {f * g}) + ({c}) * ({d * h} - {e * g})")
print(f"= {a} * ({e * i - f * h}) - ({b}) * ({d * i - f * g}) + ({c}) * ({d * h - e * g})")
print(f"= {a * (e * i - f * h)} - ({b * (d * i - f * g)}) + ({c * (d * h - e * g)})")
print(f"= {a * (e * i - f * h) + (-1) * b * (d * i - f * g) + c * (d * h - e * g)}")
print(f"= {determinant}")
