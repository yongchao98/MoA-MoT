import numpy as np

# The task is to compute the determinant of the given 3x3 matrix A.
# The matrix is the adjacency matrix of the Markov quiver for m_{1/2} = 13.
#
# A = [[ 2, -1, -1],
#      [-1,  2, -1],
#      [-1, -1,  2]]
#
# We will use the formula for the determinant of a 3x3 matrix:
# det([[a, b, c], [d, e, f], [g, h, i]]) = a(ei - fh) - b(di - fg) + c(dh - eg)

# Step 1: Define the matrix using numpy.
A = np.array([
    [2, -1, -1],
    [-1, 2, -1],
    [-1, -1, 2]
])

# Extract the elements to display the formula.
a, b, c = A[0, 0], A[0, 1], A[0, 2]
d, e, f = A[1, 0], A[1, 1], A[1, 2]
g, h, i = A[2, 0], A[2, 1], A[2, 2]

# Step 2: Calculate the terms of the determinant formula.
# First term: a(ei - fh)
term1_val = a * (e * i - f * h)
# Second term: -b(di - fg)
term2_val = -b * (d * i - f * g)
# Third term: c(dh - eg)
term3_val = c * (d * h - e * g)

# Calculate the final determinant.
determinant = term1_val + term2_val + term3_val

# Step 3: Print the step-by-step calculation.
print(f"The given matrix A is:")
print(A)
print("\nThe formula for the determinant is: det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)")

print("\nSubstituting the values from matrix A into the formula:")
print(f"det(A) = ({a}) * (({e})*({i}) - ({f})*({h})) - ({b}) * (({d})*({i}) - ({f})*({g})) + ({c}) * (({d})*({h}) - ({e})*({g}))")

print("\nCalculating the sub-expressions (the cofactors):")
print(f"det(A) = ({a}) * ({e*i - f*h}) - ({b}) * ({d*i - f*g}) + ({c}) * ({d*h - e*g})")

print("\nPerforming the multiplications:")
print(f"det(A) = ({term1_val}) + ({term2_val}) + ({term3_val})")

print("\nSumming the terms for the final result:")
print(f"det(A) = {determinant}")
