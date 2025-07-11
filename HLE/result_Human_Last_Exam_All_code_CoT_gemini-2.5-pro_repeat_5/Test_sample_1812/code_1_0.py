import numpy as np

# Define the matrix from the problem statement
A = np.array([
    [2, -1, -1],
    [-1, 2, -1],
    [-1, -1, 2]
])

# Calculate the determinant
determinant = np.linalg.det(A)

# Print the matrix
print("The matrix A is:")
print(A)
print("\nCalculating the determinant using cofactor expansion along the first row:")

# Extract elements for clarity in the printout
a11, a12, a13 = A[0, 0], A[0, 1], A[0, 2]
a21, a22, a23 = A[1, 0], A[1, 1], A[1, 2]
a31, a32, a33 = A[2, 0], A[2, 1], A[2, 2]

# Calculate intermediate terms for the equation
term1 = a11 * (a22 * a33 - a23 * a32)
term2 = -a12 * (a21 * a33 - a23 * a31)
term3 = a13 * (a21 * a32 - a22 * a31)

# Print the full equation with numbers
print(f"det(A) = ({a11}) * (({a22})*({a33}) - ({a23})*({a32})) + ({-a12}) * (({a21})*({a33}) - ({a23})*({a31})) + ({a13}) * (({a21})*({a32}) - ({a22})*({a31}))")
print(f"       = ({a11}) * ({a22 * a33 - a23 * a32}) + ({-a12}) * ({a21 * a33 - a23 * a31}) + ({a13}) * ({a21 * a32 - a22 * a31})")
print(f"       = {term1} + {term2} + {term3}")
print(f"       = {term1 + term2 + term3}")

# Print the final result
print(f"\nThe determinant of the adjacency matrix is: {determinant:.0f}")