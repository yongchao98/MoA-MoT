import numpy as np

def solve_vest_commuting(T_matrices, v, k):
    """
    Solves an instance of the VEST problem for commuting matrices.

    This function calculates the sum of products of all k-sized subsets of
    the given transformation matrices, and applies the result to vector v.
    This is FPT in k.

    Args:
        T_matrices (list of np.ndarray): A list of commuting rational matrices.
        v (np.ndarray): The initial rational vector.
        k (int): The parameter k, the size of the subsets of transformations.

    Returns:
        np.ndarray: The final vector after evaluation.
    """
    if not T_matrices:
        return np.zeros_like(v)
    
    m = len(T_matrices)
    n = T_matrices[0].shape[0]

    # M[i][j] will store the sum of products of j matrices from the first i matrices
    # M[i][j] = sum_{I subset {1..i}, |I|=j} prod_{l in I} T_l
    # We use a space-optimized DP, only keeping track of the previous column (j-1) and current column (j)
    # of the M_i,j table.
    
    # Initialize DP table for j=0. M[i][0] is always the identity matrix.
    M_prev_j = [np.identity(n) for _ in range(m + 1)]
    M_curr_j = [np.zeros((n, n)) for _ in range(m + 1)]

    for j in range(1, k + 1):
        for i in range(1, m + 1):
            # Recurrence: M[i][j] = M[i-1][j] + M[i-1][j-1] * T_i
            # In our space-optimized version:
            # M_curr_j[i] = M_curr_j[i-1] + M_prev_j[i-1] @ T_matrices[i-1]
            T_i = T_matrices[i-1]
            M_curr_j[i] = M_curr_j[i-1] + M_prev_j[i-1] @ T_i
        
        # Current column becomes previous column for the next iteration
        M_prev_j = M_curr_j
        M_curr_j = [np.zeros((n, n)) for _ in range(m + 1)]


    # The final sum of matrix products is M[m][k]
    sum_of_products_matrix = M_prev_j[m]
    
    # The problem asks for S * (sum_matrix) * v. We assume S=I.
    final_vector = sum_of_products_matrix @ v

    print("Demonstration of the FPT algorithm for VEST with commuting matrices (Part a):")
    print(f"Parameter k = {k}")
    print("\nInput matrices T_i:")
    for i, T in enumerate(T_matrices):
        print(f"T_{i+1}:\n{T}")
    
    print(f"\nInput vector v:\n{v}")
    
    print(f"\nResult of the sum of products of all {k}-subsets of matrices (M_m,k):")
    print(sum_of_products_matrix)
    
    print("\nFinal equation: M_m,k * v = Final Vector")
    # The prompt asks to "output each number in the final equation"
    # We will print the equation in a representative format for the first component
    print("Example for the first component of the final vector:")
    eq_str = " * ".join([f"{val:.2f}" for val in sum_of_products_matrix[0,:]]) + " * " + " + ".join([f"{val:.2f}" for val in v])
    # This is not quite right. Let's show the matrix-vector multiplication.
    for r in range(n):
        row_str = " ".join([f"{sum_of_products_matrix[r, c]:>5.2f}" for c in range(n)])
        v_val = v[r]
        res_val = final_vector[r]
        if r == n // 2:
            print(f"| {row_str} | * | {v_val:>5.2f} | = | {res_val:>5.2f} |")
        else:
            print(f"| {row_str} |   | {v_val:>5.2f} |   | {res_val:>5.2f} |")

    print(f"\nFinal evaluated vector:\n{final_vector}")


if __name__ == '__main__':
    # Example instance from the thought process
    # n = 2, m = 3, k = 2
    T1 = np.array([[1, 1], [0, 1]], dtype=float)
    T2 = np.array([[1, 2], [0, 1]], dtype=float)
    T3 = np.array([[2, 0], [0, 2]], dtype=float)
    
    v_vec = np.array([1, 1], dtype=float)
    k_param = 2
    
    # The matrices T1, T2, T3 commute with each other.
    # T1*T2 = [[1,3],[0,1]], T2*T1 = [[1,3],[0,1]]
    # T1*T3 = [[2,2],[0,2]], T3*T1 = [[2,2],[0,2]]
    # T2*T3 = [[2,4],[0,2]], T3*T2 = [[2,4],[0,2]]
    
    solve_vest_commuting([T1, T2, T3], v_vec, k_param)
