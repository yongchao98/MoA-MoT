import numpy as np

# The problem provides the adjacency matrix A.
A = np.array([
    [2, -1, -1],
    [-1, 2, -1],
    [-1, -1, 2]
])

# We need to compute the determinant of this 3x3 matrix.
# Let's assign the elements to variables for clarity, following the standard formula:
# | a b c |
# | d e f |
# | g h i |
a, b, c = A[0, 0], A[0, 1], A[0, 2]
d, e, f = A[1, 0], A[1, 1], A[1, 2]
g, h, i = A[2, 0], A[2, 1], A[2, 2]

# The determinant formula is: a(ei - fh) - b(di - fg) + c(dh - eg)
# Let's perform the calculation.
determinant = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)

# As requested, we will print the equation with all the numbers from the matrix.
print(f"To compute the determinant of the matrix:\n{A}\n")
print("We use the formula: det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)")
print("\nSubstituting the numbers into the equation gives:")
# This line prints the full equation with every number from the matrix.
print(f"det(A) = ({a}) * (({e})*({i}) - ({f})*({h})) - ({b}) * (({d})*({i}) - ({f})*({g})) + ({c}) * (({d})*({h}) - ({e})*({g}))")

# Show the step-by-step simplification of the equation.
term1 = e * i - f * h
term2 = d * i - f * g
term3 = d * h - e * g
print(f"det(A) = ({a})*({term1}) - ({b})*({term2}) + ({c})*({term3})")

term1_val = a * term1
term2_val = b * term2
term3_val = c * term3
print(f"det(A) = {term1_val} - ({term2_val}) + ({term3_val})")

print(f"det(A) = {term1_val + (-1 * term2_val) + term3_val}")
print(f"\nThe final determinant is: {determinant}")