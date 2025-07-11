import numpy as np

# Based on the analysis, the smallest n for f(n) > 10 is 11.
n = 11

# The Weyr matrix W_n for a matrix with a single eigenvalue and a single
# Weyr/Jordan block is a bidiagonal matrix of size n x n.
# It has 1s on the main diagonal and 1s on the first superdiagonal.
W_n = np.diag(np.ones(n)) + np.diag(np.ones(n - 1), k=1)

# Calculate the infinity norm of W_n, which is the maximum absolute row sum.
norm_inf = np.linalg.norm(W_n, ord=np.inf)

# Calculate the final result.
result = n * norm_inf

# Print the final equation as requested.
# The norm is a float, so we format it as an integer.
print(f"The value of n is: {n}")
print(f"The matrix W_n for n={n} is:")
print(W_n)
print(f"\nThe infinity norm of W_{n} is: {int(norm_inf)}")
print(f"\nThe final result is n * ||W_n||_inf:")
print(f"{n} * {int(norm_inf)} = {int(result)}")
