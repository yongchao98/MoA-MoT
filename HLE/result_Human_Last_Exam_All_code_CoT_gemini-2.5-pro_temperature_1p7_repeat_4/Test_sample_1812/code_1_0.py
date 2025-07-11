import numpy as np

# Define the adjacency matrix A
A = np.array([
    [2, -1, -1],
    [-1, 2, -1],
    [-1, -1, 2]
])

print("The given matrix A is:")
print(A)
print("\nTo compute the determinant, we use the cofactor expansion formula:")
print("det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)")

# Extract the elements of the matrix
a, b, c = A[0, 0], A[0, 1], A[0, 2]
d, e, f = A[1, 0], A[1, 1], A[1, 2]
g, h, i = A[2, 0], A[2, 1], A[2, 2]

# Print the formula with the numbers from the matrix A
print("\nPlugging in the numbers from matrix A:")
# First term
part1 = f"({a}) * (({e})*({i}) - ({f})*({h}))"
# Second term (note the minus in the formula)
part2 = f"({b}) * (({d})*({i}) - ({f})*({g}))"
# Third term
part3 = f"({c}) * (({d})*({h}) - ({e})*({g}))"

print(f"det(A) = {part1} - {part2} + {part3}")

# Calculate the intermediate values
val1 = (e * i) - (f * h)
val2 = (d * i) - (f * g)
val3 = (d * h) - (e * g)

print(f"       = ({a})*({val1}) - ({b})*({val2}) + ({c})*({val3})")

# Calculate the final result
determinant = a * val1 - b * val2 + c * val3

print(f"       = {a * val1} - {b * val2} + {c * val3}")
print(f"       = {determinant}")

# You can also use numpy's built-in function to verify
# determinant_numpy = np.linalg.det(A)

print(f"\nThe final computed determinant is: {determinant}")
