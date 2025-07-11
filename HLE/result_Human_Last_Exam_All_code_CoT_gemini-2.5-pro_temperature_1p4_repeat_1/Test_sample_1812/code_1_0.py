import numpy as np

# The matrix A is defined as given in the problem.
A = np.array([
    [2, -1, -1],
    [-1, 2, -1],
    [-1, -1, 2]
])

# Extract the elements of the matrix to display them in the equation.
a, b, c = A[0, 0], A[0, 1], A[0, 2]
d, e, f = A[1, 0], A[1, 1], A[1, 2]
g, h, i = A[2, 0], A[2, 1], A[2, 2]

# The determinant formula is det(A) = a(ei - fh) - b(di - fg) + c(dh - eg).
# We print the formula with the specific numbers from matrix A.
print("Calculating the determinant using the formula: a(ei - fh) - b(di - fg) + c(dh - eg)")
print("Substituting the values from matrix A:")
# Print the equation with each number explicitly shown.
print(f"det(A) = ({a}) * (({e})*({i}) - ({f})*({h})) - ({b}) * (({d})*({i}) - ({f})*({g})) + ({c}) * (({d})*({h}) - ({e})*({g}))")

# Calculate the intermediate terms to show the process.
term1 = e*i - f*h
term2 = d*i - f*g
term3 = d*h - e*g
print(f"det(A) = ({a}) * ({term1}) - ({b}) * ({term2}) + ({c}) * ({term3})")

# Calculate the final value of each main term.
final_term1 = a * term1
final_term2 = -b * term2
final_term3 = c * term3
print(f"det(A) = ({final_term1}) + ({final_term2}) + ({final_term3})")

# Compute the final determinant of the matrix using numpy, which is efficient and accurate.
# The result is rounded to handle potential floating point inaccuracies.
determinant = round(np.linalg.det(A))

print("\nFinal Result:")
print(f"The determinant of the adjacency matrix is: {determinant}")