import numpy as np

# The problem is to compute the determinant of the adjacency matrix of the Markov quiver.
# The context "m_1/2 = 13" is specific to the theory of Markov quivers but is not
# needed for the direct calculation of the determinant of the provided matrix.

# The given matrix A is:
A = np.array([
    [2, -1, -1],
    [-1, 2, -1],
    [-1, -1, 2]
])

# To compute the determinant of a 3x3 matrix, we use the formula:
# det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)
# where the matrix elements are:
# [[a, b, c],
#  [d, e, f],
#  [g, h, i]]

# Extract the elements from matrix A
a, b, c = A[0, 0], A[0, 1], A[0, 2]
d, e, f = A[1, 0], A[1, 1], A[1, 2]
g, h, i = A[2, 0], A[2, 1], A[2, 2]

# Calculate the individual terms
term1 = e * i - f * h
term2 = d * i - f * g
term3 = d * h - e * g

# Calculate the final determinant
determinant = a * term1 - b * term2 + c * term3

# As requested, here is the equation with all the numbers
print(f"Computing the determinant of matrix A:\n{A}\n")
print("Using the formula: det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)\n")
print("Step-by-step calculation:")

# The final equation with each number explicitly shown
print(f"det(A) = ({a}) * (({e})*({i}) - ({f})*({h})) - ({b}) * (({d})*({i}) - ({f})*({g})) + ({c}) * (({d})*({h}) - ({e})*({g}))")

# Intermediate calculation steps
print(f"       = {a} * ({e*i} - {f*h}) - ({b}) * ({d*i} - {f*g}) + ({c}) * ({d*h - e*g})")
print(f"       = {a} * ({term1}) - ({b}) * ({term2}) + {c} * ({term3})")
print(f"       = {a * term1} - ({b * term2}) + ({c * term3})")
print(f"       = {a * term1 + (-1)*(b*term2) + c*term3}")
print(f"       = {int(determinant)}")

print(f"\nThe determinant of the adjacency matrix is: {int(determinant)}")