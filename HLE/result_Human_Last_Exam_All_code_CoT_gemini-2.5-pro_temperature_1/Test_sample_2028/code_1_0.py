import numpy as np

def solve_vest_commuting(u, v, T_matrices, k):
    """
    Solves the VEST problem for commuting matrices using a dynamic programming FPT algorithm.

    Args:
        u (np.ndarray): The final evaluation vector (size n).
        v (np.ndarray): The initial vector (size n).
        T_matrices (list of np.ndarray): A list of m commuting n x n matrices.
        k (int): The parameter k.

    Returns:
        float: The result of the VEST computation.
    """
    m = len(T_matrices)
    if m == 0:
        return 0.0
    n = T_matrices[0].shape[0]

    # DP table storing vectors. dp_vectors[l] corresponds to E_l * v
    # We only need to store the previous column to compute the next one.
    # dp_vectors_prev[l] will store E_l^{(j-1)} * v
    # dp_vectors_curr[l] will store E_l^{(j)} * v
    dp_vectors_prev = [np.zeros(n) for _ in range(k + 1)]
    dp_vectors_prev[0] = v.astype(np.float64) # E_0 is the identity matrix, so E_0 * v = v

    # Iterate through each matrix T_j
    for j in range(m):
        T_j = T_matrices[j]
        dp_vectors_curr = [np.zeros(n) for _ in range(k + 1)]
        # E_0 is always I, so the vector is always v
        dp_vectors_curr[0] = v.astype(np.float64) 
        
        # Apply the recurrence E_l^{(j)} = E_l^{(j-1)} + T_j * E_{l-1}^{(j-1)}
        # This translates to v_{l,j} = v_{l,j-1} + T_j * v_{l-1,j-1}
        for l in range(1, k + 1):
            v_l_jm1 = dp_vectors_prev[l]
            v_lm1_jm1 = dp_vectors_prev[l-1]
            # The calculation step
            dp_vectors_curr[l] = v_l_jm1 + T_j @ v_lm1_jm1
        
        dp_vectors_prev = dp_vectors_curr

    # The final vector we need is dp_vectors_prev[k], which is E_k^{(m)} * v
    final_vector = dp_vectors_prev[k]
    
    # The result is the inner product u^T * final_vector
    result = u.dot(final_vector)
    
    # "output each number in the final equation!" is interpreted as printing the final sum.
    print(f"The final sum is: {result}")
    
    return result

if __name__ == '__main__':
    # Define a sample problem instance for part (a)
    n = 3  # Dimension of the vector space
    m = 4  # Number of matrices
    k = 2  # Parameter k

    # The vector u
    u = np.array([1, 1, 1])
    # The vector v
    v = np.array([1, 2, 3])

    # A list of m=4 commuting matrices T_i. Diagonal matrices always commute.
    T_matrices = [
        np.diag([1.0, 0.0, 2.0]),
        np.diag([3.0, 1.0, 0.0]),
        np.diag([1.0, 1.0, 1.0]),
        np.diag([2.0, 3.0, 2.0])
    ]

    print(f"Solving VEST for a sample instance:")
    print(f"n={n}, m={m}, k={k}")
    
    solve_vest_commuting(u, v, T_matrices, k)
