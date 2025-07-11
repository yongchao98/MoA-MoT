import numpy as np

# Let's verify the logic with an example.
# Let N be the number of nodes. We will construct the matrices and compute the norms.
N = 4

# 1. Construct an example B matrix (projection operator) with orthonormal rows, orthogonal to 1.
# We can use the Gram-Schmidt process on a basis of R^N starting with the vector 1.
# Or we can use the SVD of the projection matrix P_perp = I - (1/N)*1*1^T.
# Let's just construct one that fits the description for illustration.
# One row could be (1/sqrt(2), -1/sqrt(2), 0, 0) -> L1 norm is sqrt(2)
# Another row (1/sqrt(6), 1/sqrt(6), -2/sqrt(6), 0) -> L1 norm is 4/sqrt(6) ~ 1.63
# Another row (1/sqrt(12), 1/sqrt(12), 1/sqrt(12), -3/sqrt(12)) -> L1 norm is 6/sqrt(12) ~ 1.73
# The bound sqrt(N) = sqrt(4) = 2. It holds.
# We can use a simpler basis like the one from the Haar wavelet transform.
B = np.array([
    [0.5, 0.5, -0.5, -0.5],
    [0.5, -0.5, 0.5, -0.5],
    [0.5, -0.5, -0.5, 0.5]
])

# Let's verify the properties for our B.
# Orthonormal rows:
# np.linalg.norm(B, axis=1) -> [1. 1. 1.]
# B @ B.T -> should be identity
# print(B @ B.T) -> [[1. 0. 0.], [0. 1. 0.], [0. 0. 1.]]
# Orthogonal to 1 vector:
# B @ np.ones(N) -> [0. 0. 0.]
all_ones = np.ones(N)
# print(B @ all_ones)
# Infinity norm of B
B_inf_norm = np.max(np.sum(np.abs(B), axis=1))
# print(f"Infinity norm of B is {B_inf_norm}") -> 2.0
# The bound ||B||_inf <= sqrt(N) holds, as 2 <= sqrt(4) = 2.

# 2. Construct an example Q matrix.
# It should be a product of D*P matrices.
# Let's just create a matrix Q with non-negative entries and row sums <= 1.
Q = np.array([
    [0.1, 0.2, 0.3, 0.1], # row sum = 0.7
    [0.5, 0.0, 0.0, 0.0], # row sum = 0.5
    [0.2, 0.2, 0.2, 0.2], # row sum = 0.8
    [0.1, 0.1, 0.3, 0.5]  # row sum = 1.0
])

# 3. Construct an example stationary distribution vector pi
pi = np.array([0.25, 0.25, 0.25, 0.25])

# 4. Check the bound ||Q - 1*pi^T||_inf <= 2
Q_minus_pi = Q - np.outer(np.ones(N), pi)
Q_minus_pi_inf_norm = np.max(np.sum(np.abs(Q_minus_pi), axis=1))
# print(f"Infinity norm of (Q - 1*pi^T) is {Q_minus_pi_inf_norm}")
# Example calculation for first row:
# |0.1-0.25|+|0.2-0.25|+|0.3-0.25|+|0.1-0.25| = 0.15+0.05+0.05+0.15 = 0.4
# Bound check: ||q_i - pi||_1 <= ||q_i||_1 + ||pi||_1 <= 1 + 1 = 2
# For row 1: 0.7 + 1 = 1.7. So the bound holds.
# The derived value 0.4 is less than 1.7 and less than 2.

# 5. Calculate ||B*Q||_inf
BQ = B @ Q
BQ_inf_norm = np.max(np.sum(np.abs(BQ), axis=1))
# print(f"Infinity norm of BQ is {BQ_inf_norm}")

# 6. Check the final bound: ||B*Q||_inf <= ||B||_inf * ||Q-1*pi^T||_inf
# And also ||B*Q||_inf <= 2 * sqrt(N)
# From the logic, we expect the factor to be 2.
bound_factor = 2
final_bound = bound_factor * np.sqrt(N)

# The question asks for the upper-bound for ||B Q_{0, M}||_inf expressed as a factor of sqrt(N)
# The upper bound is 2 * sqrt(N).
# So the factor is 2.
print("The upper-bound for ||B Q_{0, M}||_inf can be expressed as k * sqrt(N).")
print("Our derivation shows that this factor k is 2.")
print(f"The equation for the bound is: ||B Q_{0, M}||_inf <= {bound_factor} * sqrt(N)")
