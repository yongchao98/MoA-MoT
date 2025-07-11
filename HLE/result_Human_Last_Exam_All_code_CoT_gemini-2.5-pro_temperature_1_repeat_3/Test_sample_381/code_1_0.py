import numpy as np

def explain_bound_calculation():
    """
    This function explains the step-by-step derivation of the upper bound for ||B Q_0,M||_inf.
    It prints the components of the final inequality.
    """

    # Step 1: Bound for ||B||_inf
    # The matrix B is the projection matrix B = I - (1/N) * 1 * 1^T.
    # The infinity norm of a matrix is the maximum absolute row sum.
    # For any row i of B, the diagonal element is (1 - 1/N) and the N-1 off-diagonal elements are -1/N.
    # The absolute row sum is |1 - 1/N| + (N-1) * |-1/N| = (1 - 1/N) + (N-1)/N = 2*(N-1)/N.
    # As N gets large, this approaches 2. So, ||B||_inf <= 2.
    bound_B = 2

    # Step 2: Bound for ||Q_0,M||_inf
    # Q_0,M is a product of matrices D_t * P_t.
    # D_t is a diagonal matrix with elements in [0, 1].
    # P_t is a row-stochastic matrix (non-negative entries, rows sum to 1).
    # The product M_t = D_t * P_t has non-negative entries.
    # The i-th row sum of M_t is d_i * sum(P_ij) = d_i * 1 = d_i, which is <= 1.
    # A matrix with non-negative entries and row sums <= 1 is sub-stochastic.
    # The product of sub-stochastic matrices is sub-stochastic.
    # Therefore, Q_0,M is sub-stochastic, and its row sums are <= 1.
    # The infinity norm of Q_0,M is its maximum absolute row sum, which is <= 1.
    bound_Q = 1

    # Step 3: Combine the bounds
    # Using the sub-multiplicative property: ||B * Q_0,M||_inf <= ||B||_inf * ||Q_0,M||_inf
    final_bound = bound_B * bound_Q

    print("The upper bound is derived from the sub-multiplicative property of matrix norms:")
    print("||B * Q_0,M||_inf <= ||B||_inf * ||Q_0,M||_inf")
    print(f"The bound for ||B||_inf is {bound_B}")
    print(f"The bound for ||Q_0,M||_inf is {bound_Q}")
    print(f"Thus, the final upper bound is {bound_B} * {bound_Q} = {final_bound}")
    print(f"The final answer is {final_bound}")

explain_bound_calculation()