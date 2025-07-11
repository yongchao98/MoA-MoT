import numpy as np

# This script demonstrates the calculation of the matrix E and its rank for a specific example.
# According to the theory, the greatest possible rank for E is 1.

# Step 1: Define an example matrix A, vector b, and a nonzero vector x.
# We choose them such that the residual r = b - Ax is a non-zero vector.
A = np.array([[1, 2, 3],
              [4, 5, 6]])
b = np.array([10, 11])
x = np.array([1, 1, 1])

print("Given an example problem with:")
print("A =\n", A)
print("b =", b)
print("x =", x)
print("-" * 30)

# Step 2: Calculate the residual vector r = b - Ax.
# The constraint on E is Ex = r.
r = b - A @ x
print("The residual r = b - Ax is calculated as:")
print("A @ x =", A @ x)
print("r = b - (A @ x) =", r)
print("-" * 30)

# Step 3: Construct the matrix E with the minimum Frobenius norm.
# The formula is E = (r * x^T) / (x^T * x).
# We will print the numbers and matrices involved in this equation.
outer_product = np.outer(r, x)
x_norm_squared = np.dot(x, x)
E = outer_product / x_norm_squared

print("The equation for E is: E = (r * x^T) / (x^T * x)\n")
print(f"The outer product (r * x^T) is:\n{outer_product}\n")
print(f"The squared norm (x^T * x) is:\n{x_norm_squared}\n")
print(f"The resulting matrix E is:\n{E}\n")
print("-" * 30)

# Step 4: Calculate the rank of the resulting matrix E.
# As r and x are non-zero, the rank must be 1.
rank_E = np.linalg.matrix_rank(E)

print(f"The rank of matrix E is: {rank_E}")
print("\nThis demonstrates that a rank of 1 is possible, which is the greatest possible rank.")
