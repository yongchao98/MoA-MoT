import numpy as np

# Set the dimensions for our demonstration
m = 4  # Number of rows for A and b
n = 3  # Number of columns for A, and rows for x

# For reproducibility, we fix the random seed
np.random.seed(0)

# 1. Define a random matrix A, and non-zero vectors b and x.
#    This setup represents a general case.
A = np.random.randn(m, n)
x = np.random.randn(n, 1)
b = np.random.randn(m, 1)

print("A =")
print(A)
print("\nx =")
print(x)
print("\nb =")
print(b)

# 2. Calculate the residual vector r = b - Ax.
#    With random inputs, r will be non-zero with near certainty.
r = b - A @ x
print("\nThe residual r = b - Ax =")
print(r)

# 3. Calculate the matrix E with the minimum Frobenius norm using the derived formula:
#    E = (r * x^T) / (x^T * x)
xtx = x.T @ x  # This is a scalar value, ||x||_2^2
E = (r @ x.T) / xtx

print("\nThe minimum-norm matrix E = (r * x^T) / (x^T * x) =")
print(E)

# 4. Calculate the rank of the resulting matrix E.
rank_E = np.linalg.matrix_rank(E)

# 5. Output the final computed rank.
#    This demonstrates the greatest possible rank derived in the explanation.
print(f"\nThe rank of the matrix E is {rank_E}.")
