import numpy as np

# The given adjacency matrix A
A = np.array([
    [2, -1, -1],
    [-1, 2, -1],
    [-1, -1, 2]
])

# Extracting the elements for clarity in the formula
a, b, c = A[0, 0], A[0, 1], A[0, 2]
d, e, f = A[1, 0], A[1, 1], A[1, 2]
g, h, i = A[2, 0], A[2, 1], A[2, 2]

# Calculate the determinant using numpy.linalg.det
# The result might be a very small float due to precision, so we round it.
determinant = round(np.linalg.det(A))

# Print the step-by-step calculation as an equation
print("The determinant is calculated using the formula: a(ei - fh) - b(di - fg) + c(dh - eg)")
print("Plugging in the values from matrix A:")
print(f"det(A) = ({a}) * (({e})*({i}) - ({f})*({h})) - ({b}) * (({d})*({i}) - ({f})*({g})) + ({c}) * (({d})*({h}) - ({e})*({g}))")

# Calculate intermediate terms
term1 = a * (e*i - f*h)
term2 = -b * (d*i - f*g)
term3 = c * (d*h - e*g)

print(f"det(A) = {term1} + {term2} + {term3}")

# Final result
print(f"det(A) = {determinant}")

print("\nThe final determinant of the adjacency matrix is:")
print(determinant)