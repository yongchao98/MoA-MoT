import numpy as np

# The problem is to compute the determinant of the matrix A.
# The context about the Markov quiver and m_1/2 = 13 leads to this specific matrix,
# but the calculation itself is a standard determinant problem.
A = np.array([
    [2, -1, -1],
    [-1, 2, -1],
    [-1, -1, 2]
])

# For a 3x3 matrix [[a, b, c], [d, e, f], [g, h, i]], the determinant is:
# det = a(ei - fh) - b(di - fg) + c(dh - eg)

# Extracting the elements for clarity
a, b, c = A[0, 0], A[0, 1], A[0, 2]
d, e, f = A[1, 0], A[1, 1], A[1, 2]
g, h, i = A[2, 0], A[2, 1], A[2, 2]

# Calculate the terms
term1 = a * (e * i - f * h)
term2 = -b * (d * i - f * g)
term3 = c * (d * h - e * g)
determinant = term1 + term2 + term3

# Print the step-by-step calculation as requested
print("The matrix A is:")
print(A)
print("\nThe determinant is calculated using the formula: det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)")

print("\nSubstituting the values from matrix A into the formula:")
print(f"det(A) = ({a}) * (({e})*({i}) - ({f})*({h})) - ({b}) * (({d})*({i}) - ({f})*({g})) + ({c}) * (({d})*({h}) - ({e})*({g}))")

print("\nCalculating the intermediate values:")
print(f"det(A) = ({a}) * ({e*i - f*h}) - ({b}) * ({d*i - f*g}) + ({c}) * ({d*h - e*g})")
print(f"det(A) = ({term1}) + ({term2}) + ({term3})")

print("\nFinal Result:")
print(f"det(A) = {determinant}")