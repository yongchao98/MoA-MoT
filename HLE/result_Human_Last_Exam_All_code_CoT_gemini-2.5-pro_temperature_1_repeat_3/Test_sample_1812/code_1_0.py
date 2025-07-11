import numpy as np

# The matrix is given as:
# A = [[ 2, -1, -1],
#      [-1,  2, -1],
#      [-1, -1,  2]]
# The information m_{1/2} = 13 is not needed for this calculation.

A = np.array([
    [2, -1, -1],
    [-1, 2, -1],
    [-1, -1, 2]
])

# Extracting the elements for clarity in the formula
a, b, c = A[0, 0], A[0, 1], A[0, 2]
d, e, f = A[1, 0], A[1, 1], A[1, 2]
g, h, i = A[2, 0], A[2, 1], A[2, 2]

# Calculate the determinant using the formula:
# det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)
term1 = a * (e * i - f * h)
term2 = b * (d * i - f * g)
term3 = c * (d * h - e * g)
determinant = term1 - term2 + term3

print(f"Given the matrix A:\n{A}\n")
print("The determinant is calculated using the formula: det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)")
print("Substituting the values from matrix A:")

# Printing each number in the final equation
print(f"det(A) = {a}(({e})*({i}) - ({f})*({h})) - ({b})(({d})*({i}) - ({f})*({g})) + ({c})(({d})*({h}) - ({e})*({g}))")
print(f"       = {a}({e*i} - {f*h}) - ({b})({d*i} - {f*g}) + ({c})({d*h} - {e*g})")
print(f"       = {a}({e*i - f*h}) - ({b})({d*i - f*g}) + ({c})({d*h - e*g})")
print(f"       = {term1} - ({term2}) + ({term3})")
print(f"       = {term1 + term2 + term3}") # Note: term2 is negative, so term1 - term2 becomes term1 + (-term2)
print(f"       = {determinant}")

# Verifying with numpy's built-in function
np_determinant = np.linalg.det(A)
print(f"\nDeterminant calculated by NumPy: {np_determinant}")
print(f"The final result is {int(round(np_determinant))}.")

<<<0>>>