import numpy as np

# The problem asks to compute the determinant of the adjacency matrix of the Markov quiver.
# The matrix A is explicitly provided. The value m_1/2 = 13 is not needed for this calculation.

# Define the matrix A
A = np.array([
    [2, -1, -1],
    [-1, 2, -1],
    [-1, -1, 2]
])

# Extract elements for clarity in the formula representation
a, b, c = A[0, 0], A[0, 1], A[0, 2]
d, e, f = A[1, 0], A[1, 1], A[1, 2]
g, h, i = A[2, 0], A[2, 1], A[2, 2]

# --- Calculation steps ---

# 1. Calculate the values of the 2x2 sub-determinants (minors)
sub_det1 = e * i - f * h
sub_det2 = d * i - f * g
sub_det3 = d * h - e * g

# 2. Calculate each term of the determinant formula
term1 = a * sub_det1
term2 = -b * sub_det2
term3 = c * sub_det3

# 3. Sum the terms to get the final determinant
determinant = term1 + term2 + term3

# --- Print the process ---
print("Given the matrix A:")
print(A)
print("\nWe compute the determinant using the formula for a 3x3 matrix:")
print("det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)\n")

# Print the equation with all the numbers from the matrix
print("Step 1: Plugging the numbers into the formula:")
print(f"det(A) = ({a}) * (({e})*({i}) - ({f})*({h})) - ({b}) * (({d})*({i}) - ({f})*({g})) + ({c}) * (({d})*({h}) - ({e})*({g}))")

# Print the equation with the inner parts calculated
print("\nStep 2: Evaluating the expressions inside the parentheses:")
print(f"det(A) = ({a}) * ({sub_det1}) - ({b}) * ({sub_det2}) + ({c}) * ({sub_det3})")

# Print the equation with the multiplications done
print("\nStep 3: Performing the multiplications:")
print(f"det(A) = ({term1}) + ({term2}) + ({term3})")

# Print the final result
print("\nStep 4: Final calculation:")
print(f"det(A) = {determinant}")