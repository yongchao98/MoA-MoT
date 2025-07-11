import numpy as np

# Define the matrix A
A = np.array([
    [2, -1, -1],
    [-1, 2, -1],
    [-1, -1, 2]
])

# Extract the elements for the equation printout
a, b, c = A[0, 0], A[0, 1], A[0, 2]
d, e, f = A[1, 0], A[1, 1], A[1, 2]
g, h, i = A[2, 0], A[2, 1], A[2, 2]

# Calculate the determinant using numpy
determinant = np.linalg.det(A)

# Print the final equation with all the numbers
print("The determinant is calculated using the formula: a(ei - fh) - b(di - fg) + c(dh - eg)")
print("Substituting the values from the matrix:")
print(f"Determinant = ({a}) * (({e})*({i}) - ({f})*({h})) - ({b}) * (({d})*({i}) - ({f})*({g})) + ({c}) * (({d})*({h}) - ({e})*({g}))")
print(f"Determinant = ({a}) * ({e*i - f*h}) + {-b} * ({d*i - f*g}) + ({c}) * ({d*h - e*g})")
print(f"Determinant = {a * (e*i - f*h)} + {-b * (d*i - f*g)} + {c * (d*h - e*g)}")
print(f"Determinant = {int(round(determinant))}")

<<<0>>>