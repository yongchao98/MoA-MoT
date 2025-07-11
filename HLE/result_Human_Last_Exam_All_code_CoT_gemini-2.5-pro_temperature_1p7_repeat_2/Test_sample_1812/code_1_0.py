import numpy as np

# Define the matrix A
A = np.array([
    [2, -1, -1],
    [-1, 2, -1],
    [-1, -1, 2]
])

# Extract the elements of the matrix for clarity in the formula
a, b, c = A[0, 0], A[0, 1], A[0, 2]
d, e, f = A[1, 0], A[1, 1], A[1, 2]
g, h, i = A[2, 0], A[2, 1], A[2, 2]

# Calculate the determinant using numpy for the final answer
determinant = np.linalg.det(A)

# Print the step-by-step calculation as requested
print("To compute the determinant of the matrix A:")
print(A)
print("\nWe use the formula: det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)")
print("\nSubstituting the values from matrix A into the formula:")
print(f"det(A) = ({a}) * (({e})*({i}) - ({f})*({h})) - ({b}) * (({d})*({i}) - ({f})*({g})) + ({c}) * (({d})*({h}) - ({e})*({g}))")
# Perform the calculations within the parentheses
ei = e * i
fh = f * h
di = d * i
fg = f * g
dh = d * h
eg = e * g
print(f"det(A) = ({a}) * ({ei} - ({fh})) - ({b}) * ({di} - ({fg})) + ({c}) * ({dh} - ({eg}))")
# Perform the subtraction within the parentheses
term1_inner = ei - fh
term2_inner = di - fg
term3_inner = dh - eg
print(f"det(A) = ({a}) * ({term1_inner}) - ({b}) * ({term2_inner}) + ({c}) * ({term3_inner})")
# Perform the final multiplications and additions/subtractions
term1 = a * term1_inner
term2 = -b * term2_inner
term3 = c * term3_inner
print(f"det(A) = {term1} + {term2} + {term3}")
print(f"det(A) = {term1 + term2 + term3}")

# Print the final answer clearly
print(f"\nThe determinant of the matrix A is {determinant}.")
<<<0>>>