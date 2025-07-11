import numpy as np

# This script demonstrates the solution by constructing an example.
# We choose A, b, and a non-zero x, then find the minimal norm E
# such that (A+E)x = b.

# 1. Define problem dimensions and variables
m = 3  # Number of rows
n = 2  # Number of columns

A = np.array([[1.0, 2.0],
              [3.0, 4.0],
              [5.0, 6.0]])
x = np.array([[10.0],
              [20.0]]) # A non-zero vector

# Choose b such that r = b - Ax is non-zero
b = np.array([[1.0],
              [1.0],
              [1.0]])

# 2. Calculate the residual vector r = b - Ax
r = b - A @ x

# 3. Calculate the minimum Frobenius norm matrix E
# The formula is E = (r * x^T) / (x^T * x)
x_transpose = x.T
r_x_transpose = r @ x_transpose
x_transpose_x = x.T @ x

E = r_x_transpose / x_transpose_x

# 4. Calculate the rank of E
rank_E = np.linalg.matrix_rank(E)

# The instruction asks to "output each number in the final equation".
# The final equation is E = (r x^T) / (x^T x). Let's print its components.

print("--- Vector r = b - A @ x ---")
# Flatten to make it easier to read
print(r.flatten())
print("\n--- Vector x^T ---")
print(x_transpose.flatten())
print("\n--- Scalar value of x^T @ x ---")
print(x_transpose_x[0, 0])
print("\n--- Matrix E = (r @ x^T) / (x^T @ x) ---")
print(E)
print("\n--- Rank of Matrix E ---")
print(int(rank_E))

print("\nAs shown, for a case where r is non-zero, the rank of E is 1.")
print("Since the solution E is always an outer product, its rank cannot exceed 1.")
print("Therefore, the greatest possible rank is 1.")
