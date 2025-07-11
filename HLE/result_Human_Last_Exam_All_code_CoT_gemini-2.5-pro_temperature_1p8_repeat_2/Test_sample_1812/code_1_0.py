import numpy as np

# Define the adjacency matrix A
A = np.array([
    [2, -1, -1],
    [-1, 2, -1],
    [-1, -1, 2]
])

print("The adjacency matrix A is:")
print(A)
print("\nTo compute the determinant, we use the formula:")
print("det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)")

# Extract the elements of the matrix
a, b, c = A[0, 0], A[0, 1], A[0, 2]
d, e, f = A[1, 0], A[1, 1], A[1, 2]
g, h, i = A[2, 0], A[2, 1], A[2, 2]

# Calculate the terms
term1_inner = e * i - f * h
term2_inner = d * i - f * g
term3_inner = d * h - e * g

term1 = a * term1_inner
term2 = b * term2_inner
term3 = c * term3_inner

# The final result
determinant = term1 - term2 + term3

# Print the step-by-step calculation with all the numbers
print("\nSubstituting the values from matrix A:")
print(f"det(A) = ({a}) * (({e})*({i}) - ({f})*({h})) - ({b}) * (({d})*({i}) - ({f})*({g})) + ({c}) * (({d})*({h}) - ({e})*({g}))")
print(f"       = ({a}) * ({e*i} - {f*h}) - ({b}) * ({d*i} - {f*g}) + ({c}) * ({d*h} - {e*g})")
print(f"       = ({a}) * ({term1_inner}) - ({b}) * ({term2_inner}) + ({c}) * ({term3_inner})")
print(f"       = ({term1}) - ({term2}) + ({term3})")
print(f"       = {term1} + {-term2} + {term3}")
print(f"       = {determinant}")

print(f"\nThe determinant of the adjacency matrix is: {determinant}")