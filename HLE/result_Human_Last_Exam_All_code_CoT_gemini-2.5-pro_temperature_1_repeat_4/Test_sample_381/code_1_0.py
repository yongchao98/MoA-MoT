import numpy as np

# The problem asks for the upper-bound for ||B * Q_{0,M}||_inf.
# From the submultiplicative property of matrix norms, we have:
# ||B * Q_{0,M}||_inf <= ||B||_inf * ||Q_{0,M}||_inf

# 1. Bounding ||Q_{0,M}||_inf
# Q_{0,M} is defined as a product of matrices: D^(M)P^(M) ... D^(0)P^(0).
# The norm of a product is less than or equal to the product of the norms:
# ||Q_{0,M}||_inf <= ||D^(M)||_inf * ||P^(M)||_inf * ... * ||D^(0)||_inf * ||P^(0)||_inf
#
# From the text:
# - P^(t) are row-stochastic matrices. The infinity norm of a row-stochastic matrix is 1.
# - D^(t) are diagonal matrices with diagonal entries between 0 and 1.
#   The infinity norm of a diagonal matrix is the maximum absolute value of its entries.
#   So, ||D^(t)||_inf <= 1.
#
# Therefore, ||Q_{0,M}||_inf <= 1 * 1 * ... * 1 * 1 = 1.

# 2. Bounding ||B||_inf
# B is the projection matrix onto the space orthogonal to span{1}.
# B can be written as B = I - (1/N) * 1 * 1^T, where I is the NxN identity matrix
# and 1 is the all-one vector of length N.
# The infinity norm of a matrix is the maximum absolute row sum.
# Let's consider the i-th row of B.
# The diagonal element is B_ii = 1 - 1/N.
# The off-diagonal elements are B_ij = -1/N for j != i.
# The absolute row sum for any row i is:
# sum_j(|B_ij|) = |1 - 1/N| + (N-1) * |-1/N|
#               = (1 - 1/N) + (N-1)/N
#               = (N - 1 + N - 1) / N
#               = (2N - 2) / N
#               = 2 - 2/N
#
# So, ||B||_inf = 2 - 2/N.
# For any N >= 1, we have 0 <= 2/N <= 2, so ||B||_inf <= 2.

# 3. Final Bound
# Combining the two bounds:
# ||B * Q_{0,M}||_inf <= ||B||_inf * ||Q_{0,M}||_inf <= (2 - 2/N) * 1 <= 2.
#
# The question asks for the upper-bound expressed as a factor of sqrt(N).
# The derived bound is a constant, 2. This bound is rigorous and derived directly
# from the provided definitions. The other information in the prompt appears to
# be for a more detailed analysis (e.g., proving convergence to 0) rather
# than finding a simple, constant upper bound. The most reasonable interpretation
# is that the question asks for this constant upper bound.
upper_bound = 2.0
# The question asks to output the number in the final equation.
# The final equation is ||B * Q_{0,M}||_inf <= 2.0
print(upper_bound)