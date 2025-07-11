import numpy as np

# The problem asks to compute the determinant of the matrix A.
A = np.array([
    [2, -1, -1],
    [-1, 2, -1],
    [-1, -1, 2]
])

# Extract elements for clarity in the explanation.
# Matrix A = [[a, b, c], [d, e, f], [g, h, i]]
a, b, c = A[0, 0], A[0, 1], A[0, 2]
d, e, f = A[1, 0], A[1, 1], A[1, 2]
g, h, i = A[2, 0], A[2, 1], A[2, 2]

# The determinant of a 3x3 matrix is calculated using the formula:
# det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)
# The script will now print out the steps of this calculation.

print("The matrix A is:")
print(A)
print("\nTo compute the determinant, we use the cofactor expansion formula:")
print("det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)\n")

print("Substituting the values from matrix A, the final equation and its calculation are:")

# Calculate the values of the cofactors
val1 = e * i - f * h
val2 = d * i - f * g
val3 = d * h - e * g

# Calculate the full terms
term1 = a * val1
term2 = -b * val2
term3 = c * val3

# Final result
final_result = term1 + term2 + term3

# Print the step-by-step calculation
print(f"det(A) = {a}(({e})*({i}) - ({f})*({h})) - ({b})(({d})*({i}) - ({f})*({g})) + ({c})(({d})*({h}) - ({e})*({g}))")
print(f"       = {a}({e*i} - {f*h}) - ({b})({d*i} - {f*g}) + ({c})({d*h} - {e*g})")
print(f"       = {a}({val1}) - ({b})({val2}) + ({c})({val3})")
print(f"       = {term1} + ({term2}) + ({term3})")
print(f"       = {final_result}")
