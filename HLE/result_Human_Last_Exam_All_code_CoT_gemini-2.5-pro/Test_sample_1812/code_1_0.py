import numpy as np

# Define the adjacency matrix A
A = np.array([
    [2, -1, -1],
    [-1, 2, -1],
    [-1, -1, 2]
])

# Calculate the determinant using numpy's linear algebra module
determinant = np.linalg.det(A)

# Extract the individual elements of the matrix to display the calculation
a, b, c = A[0, 0], A[0, 1], A[0, 2]
d, e, f = A[1, 0], A[1, 1], A[1, 2]
g, h, i = A[2, 0], A[2, 1], A[2, 2]

# Print the calculation steps
print("The given matrix is:")
print(A)
print("\nThe determinant is calculated using the formula: a(ei - fh) - b(di - fg) + c(dh - eg)")
print("\nSubstituting the values from the matrix:")
print(f"Determinant = ({a}) * (({e})*({i}) - ({f})*({h})) - ({b}) * (({d})*({i}) - ({f})*({g})) + ({c}) * (({d})*({h}) - ({e})*({g}))")

# Calculate the value of each sub-determinant
term1_val = (e * i) - (f * h)
term2_val = (d * i) - (f * g)
term3_val = (d * h) - (e * g)
print(f"Determinant = ({a})*({term1_val}) - ({b})*({term2_val}) + ({c})*({term3_val})")

# Calculate the final terms of the sum
final_term1 = a * term1_val
final_term2 = -b * term2_val
final_term3 = c * term3_val
print(f"Determinant = {final_term1} + {final_term2} + {final_term3}")

# Print the final result, rounding to handle potential floating point inaccuracies
print(f"\nThe final computed determinant is: {int(round(determinant))}")