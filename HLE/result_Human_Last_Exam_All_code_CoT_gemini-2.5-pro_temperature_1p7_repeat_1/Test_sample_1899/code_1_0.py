import numpy as np

# The matrix S is derived from the power sums of the roots of x^3 + x^2 + x - 1 = 0.
# p_k = sum(root_i^k)
# p_0 = 3
# p_1 = -1
# p_2 = -1
# p_3 = 5
# p_4 = -5
# S = [[p_0, p_1, p_2], [p_1, p_2, p_3], [p_2, p_3, p_4]]
S = np.array([
    [3., -1., -1.],
    [-1., -1., 5.],
    [-1., 5., -5.]
])

# Find the eigenvalues of S. Since S is real and symmetric, eigenvalues are real.
eigenvalues = np.linalg.eigvalsh(S)

# Sort eigenvalues to easily find min and max
eigenvalues.sort()
lambda_min = eigenvalues[0]
lambda_max = eigenvalues[2]

# The infimum and supremum are given by 2/lambda_max and 2/lambda_min respectively.
inf_value = 2 / lambda_max
sup_value = 2 / lambda_min

# Print the results
print(f"inf_n |P_n(xi)| (a_n^2+b_n^2+c_n^2) = {inf_value}")
print(f"sup_n |P_n(xi)| (a_n^2+b_n^2+c_n^2) = {sup_value}")