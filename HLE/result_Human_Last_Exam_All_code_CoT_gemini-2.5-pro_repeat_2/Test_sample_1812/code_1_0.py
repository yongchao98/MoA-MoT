import numpy as np

# Define the matrix A
A = np.array([
    [2, -1, -1],
    [-1, 2, -1],
    [-1, -1, 2]
])

# Extract the elements for the equation display
a, b, c = A[0, 0], A[0, 1], A[0, 2]
d, e, f = A[1, 0], A[1, 1], A[1, 2]
g, h, i = A[2, 0], A[2, 1], A[2, 2]

# Calculate the determinant using numpy
determinant = np.linalg.det(A)

# Print the explanation and the calculation steps
print(f"The given matrix is:\n{A}\n")
print("The determinant is calculated using the formula: a(ei - fh) - b(di - fg) + c(dh - eg)")
print("Substituting the values from the matrix:")

# Print the full equation with numbers
print(f"Determinant = ({a}) * (({e})*({i}) - ({f})*({h})) - ({b}) * (({d})*({i}) - ({f})*({g})) + ({c}) * (({d})*({h}) - ({e})*({g}))")

# Print the intermediate calculation
term1 = e * i - f * h
term2 = d * i - f * g
term3 = d * h - e * g
print(f"Determinant = ({a}) * ({term1}) - ({b}) * ({term2}) + ({c}) * ({term3})")

# Print the final calculation step
final_calc = a * term1 - b * term2 + c * term3
print(f"Determinant = {a * term1} - {b * term2} + {c * term3}")
print(f"Determinant = {final_calc}")

print("\nFinal computed determinant:")
# Use round() to handle potential floating point inaccuracies for a clean integer result
print(round(determinant))