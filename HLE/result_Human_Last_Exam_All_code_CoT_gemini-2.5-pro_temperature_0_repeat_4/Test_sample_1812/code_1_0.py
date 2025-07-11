import numpy as np

# The task is to compute the determinant of the given matrix A.
# The context about "Markov quiver for m_{1/2} = 13" is additional information
# that does not change the matrix A, which is explicitly provided.

# Step 1: Define the matrix A.
A = np.array([
    [2, -1, -1],
    [-1, 2, -1],
    [-1, -1, 2]
])

# Step 2: Extract the elements of the matrix for the formula.
# The formula for the determinant of a 3x3 matrix [[a, b, c], [d, e, f], [g, h, i]] is:
# det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)
a, b, c = A[0, 0], A[0, 1], A[0, 2]
d, e, f = A[1, 0], A[1, 1], A[1, 2]
g, h, i = A[2, 0], A[2, 1], A[2, 2]

# Step 3: Print the calculation process step-by-step.
print(f"The given matrix A is:\n{A}\n")
print("The determinant is calculated using the formula: a(ei - fh) - b(di - fg) + c(dh - eg)")
print("Substituting the values from matrix A into the formula:")

# Print the formula with the numbers substituted
print(f"det(A) = ({a}) * (({e}) * ({i}) - ({f}) * ({h})) - ({b}) * (({d}) * ({i}) - ({f}) * ({g})) + ({c}) * (({d}) * ({h}) - ({e}) * ({g}))")

# Calculate the values inside the parentheses (the cofactors)
val1 = e * i - f * h
val2 = d * i - f * g
val3 = d * h - e * g
print(f"       = {a} * ({val1}) - ({b}) * ({val2}) + ({c}) * ({val3})")

# Calculate the products
prod1 = a * val1
prod2 = -b * val2
prod3 = c * val3
print(f"       = {prod1} + ({prod2}) + ({prod3})")

# Calculate the final sum
determinant = prod1 + prod2 + prod3
print(f"       = {determinant}")

# As a verification, we can use numpy's built-in function
# Note: numpy might return a very small floating point number close to 0 due to precision.
final_determinant = round(np.linalg.det(A))
print(f"\nThe determinant of the adjacency matrix is: {final_determinant}")