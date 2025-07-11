import math

# The problem is to find the least upper bound of the product E_P * E_H * S_P * S_H over all positive integers n.
# Based on the analysis, the product can be expressed as a function of n: F(n) = 2 + 2/n.
# This function is monotonically decreasing for n > 0.
# The least upper bound (supremum) for positive integers n (n>=1) is the value at n=1.

# We will calculate the components of the product for n=1 to demonstrate the result.

# Set n, the dimension of the simplex
n = 1

# Size of the Cayley-Menger matrix C_n is N x N, where N = n + 1
N = n + 1

# --- Calculations for H ---
# Eigenvalues of C_n, and thus H, are n and -1.
lambda_max_H = n
lambda_min_H = -1
# E_H is the average eigenvalue gap of H.
E_H = (lambda_max_H - lambda_min_H) / (N - 1)

# S_H is the mean square of the singular values of H.
# This is equal to S_C, which simplifies to n.
S_H = n

# --- Calculations for P ---
# P is a symmetric orthogonal matrix of size N x N. Its eigenvalues are +1 and -1.
lambda_max_P = 1
lambda_min_P = -1
# E_P is the average eigenvalue gap of P.
E_P = (lambda_max_P - lambda_min_P) / (N - 1)

# S_P is the mean square of the singular values of P.
# For any orthogonal matrix, all singular values are 1.
S_P = 1.0

# --- Final Product Calculation ---
# The product is F(n) = E_P * E_H * S_P * S_H
product = E_P * E_H * S_P * S_H

# The LUB is the value of the product at n=1.
LUB = product

print(f"The calculation for n = {n}:")
print(f"E_P = ({lambda_max_P} - ({lambda_min_P})) / ({N} - 1) = {E_P}")
print(f"E_H = ({lambda_max_H} - ({lambda_min_H})) / ({N} - 1) = {E_H}")
print(f"S_P = {S_P}")
print(f"S_H = {S_H}")
print(f"The product E_P * E_H * S_P * S_H for n={n} is: {E_P} * {E_H} * {S_P} * {S_H} = {product}")
print(f"The function F(n) = 2 + 2/n is decreasing, so the least upper bound is the value at n=1.")
print(f"Least Upper Bound = {LUB}")
<<<4>>>