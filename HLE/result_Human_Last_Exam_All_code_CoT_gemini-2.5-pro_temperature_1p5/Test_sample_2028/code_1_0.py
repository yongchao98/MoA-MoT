import numpy as np

def solve_vest_decision_c(v, T_matrices, S, k, t):
    """
    Solves the decision version of the VEST problem under the restriction
    for question (c), demonstrating it is in FPT (actually, P).

    The problem is to decide if sum_{I, |I|=k} S * (prod_{i in I} T_i) * v == t,
    where the product is taken over indices in increasing order.

    This function uses dynamic programming to compute the sum of matrix products.
    """
    if not T_matrices:
        m = 0
    else:
        m = len(T_matrices)
    
    if k < 0 or k > m:
        # No subsets of size k, sum is zero vector.
        n = S.shape[0]
        final_sum_matrix = np.zeros((n, n))
    elif k == 0:
        # One subset of size 0 (the empty set), product is identity.
        n = S.shape[0]
        final_sum_matrix = np.identity(n)
    else:
        n = T_matrices[0].shape[0]
        # dp[i][j] will store the sum of products for subsets of size j
        # from the first i matrices {T_0, ..., T_{i-1}}.
        # Each element dp[i][j] is an n x n matrix.
        dp = [[np.zeros((n, n)) for _ in range(k + 1)] for _ in range(m + 1)]

        # Base case: The sum over subsets of size 0 is the identity matrix
        # (product over the empty set).
        for i in range(m + 1):
            dp[i][0] = np.identity(n)

        # Fill the DP table
        for i in range(1, m + 1):
            T_prev = T_matrices[i-1]
            for j in range(1, k + 1):
                # A(i, j) = A(i-1, j) + A(i-1, j-1) @ T_{i-1}
                dp[i][j] = dp[i-1][j] + dp[i-1][j-1] @ T_prev
        
        final_sum_matrix = dp[m][k]

    # Compute the final vector
    final_vector = S @ final_sum_matrix @ v
    
    print("--- VEST Calculation Result ---")
    print(f"Final computed vector:\n{final_vector}")
    print(f"Target vector t:\n{t}")

    # The final equation to check
    # We print each number by iterating through the numpy arrays.
    print("\nIs the final equation satisfied?")
    
    is_equal = np.allclose(final_vector, t)
    
    # Building the string for the equation
    final_vector_str = "[" + ", ".join([f"{x:.4f}" for x in final_vector]) + "]"
    t_vector_str = "[" + ", ".join([f"{x:.4f}" for x in t]) + "]"
    
    print(f"{final_vector_str} == {t_vector_str}  ?  ->  {is_equal}")

    return is_equal

if __name__ == '__main__':
    # Example usage for the function
    
    # Parameters for the VEST instance
    n = 2  # Dimension of the vector space
    m = 4  # Number of matrices T_i
    k = 2  # Size of the subsets I

    # Define the input vector v
    v = np.array([1.0, 2.0])

    # Define the transformation matrices T_i.
    # Each T_i has exactly one non-zero entry per row, as per restriction (c).
    T_matrices = [
        np.array([[0.0, 1.5], [2.0, 0.0]]),
        np.array([[1.0, 0.0], [0.0, 0.5]]),
        np.array([[0.0, -1.0], [-1.0, 0.0]]),
        np.array([[0.5, 0.0], [0.0, 2.0]])
    ]

    # Define the S matrix
    S = np.array([[1.0, 0.0], [0.0, 1.0]])

    # Define a target vector t
    t = np.array([-1.5, 9.5])
    
    # Run the solver
    solve_vest_decision_c(v, T_matrices, S, k, t)
