import numpy as np

# Based on the analysis, we need to find the smallest integer n > 10.
n = 11

# The matrix W_n is the Jordan block J_n(1) because the underlying Schur matrix S_n
# is non-derogatory (since its c_1 Taylor coefficient is non-zero).
# We construct this n x n matrix, W_n.
W_n = np.diag(np.ones(n)) + np.diag(np.ones(n - 1), 1)

# Calculate the infinity norm of W_n, which is the maximum absolute row sum.
# For a Jordan block J_n(1), this norm is always 2 (for n > 1).
norm_W_n = np.linalg.norm(W_n, np.inf)

# Calculate the final result by multiplying n with the norm.
result = n * norm_W_n

# As requested, print the final equation showing all the numbers involved.
print(f"{n} * {norm_W_n} = {result}")

# Finally, output the answer in the specified format.
print(f"<<<{result}>>>")