import numpy as np

def find_null_space_Z2(M):
    """
    Finds a basis for the null space of a matrix M over Z_2 using Gaussian elimination.
    """
    rows, cols = M.shape
    A = np.copy(M) % 2

    # To row echelon form
    pivot_row = 0
    pivot_cols = []
    for j in range(cols):
        if pivot_row >= rows:
            break
        i = pivot_row
        while i < rows and A[i, j] == 0:
            i += 1
        
        if i < rows:
            A[[pivot_row, i]] = A[[i, pivot_row]] # Swap rows
            pivot_cols.append(j)
            
            for k in range(rows):
                if k != pivot_row and A[k, j] == 1:
                    A[k, :] = (A[k, :] + A[pivot_row, :]) % 2
            pivot_row += 1

    # Find basis from RREF
    null_space_basis = []
    free_cols = [j for j in range(cols) if j not in pivot_cols]
    
    for free_col in free_cols:
        v = np.zeros(cols, dtype=np.int64)
        v[free_col] = 1
        for i, pivot_col in enumerate(pivot_cols):
            v[pivot_col] = A[i, free_col]
        null_space_basis.append(v)
        
    return np.array(null_space_basis)

def solve_homogeneous_binary(A, q):
    """
    Finds a non-zero vector x in {0,1}^m such that Ax = 0 (mod q),
    where q = 2^k.

    Args:
        A (numpy.ndarray): The n x m matrix.
        q (int): The modulus, must be a power of 2, q = 2^k with k > 1.

    Returns:
        numpy.ndarray: A non-zero m-dimensional vector x, or None if no solution is found.
    """
    if q <= 1 or (q & (q - 1)) != 0:
        raise ValueError("q must be a power of 2 greater than 1.")
        
    n, m = A.shape
    k = int(np.log2(q))

    # Step 1: Base case, find null space basis for A mod 2
    A_mod_2 = A % 2
    current_basis = find_null_space_Z2(A_mod_2).T
    
    if current_basis.shape[1] == 0:
        print("No non-trivial solution modulo 2 found.")
        return None

    # Lifting steps
    modulus = 2
    for j in range(1, k):
        modulus *= 2
        
        if current_basis.shape[1] == 0:
            print(f"Algorithm failed at lifting step to modulo {modulus}: no basis vectors.")
            return None
            
        # For each basis vector u_i, calculate y_i = (A @ u_i) / 2^j
        Y_cols = []
        for i in range(current_basis.shape[1]):
            u_i = current_basis[:, i]
            # Ensure calculations are done with large enough integers
            y_i = (A.astype(np.int64) @ u_i) // (modulus // 2)
            Y_cols.append(y_i)
        
        # Form matrix Y and find its null space over Z_2
        Y = np.array(Y_cols).T
        coeffs_basis = find_null_space_Z2(Y)
        
        if coeffs_basis.shape[0] == 0:
             print(f"Algorithm failed at lifting step to modulo {modulus}: no non-trivial combination found.")
             return None

        # Create new basis for the next level
        # new_basis_vectors = (current_basis @ coeffs_basis.T) % 2
        new_basis_vectors = np.dot(current_basis, coeffs_basis.T) % 2
        current_basis = new_basis_vectors

    if current_basis.shape[1] == 0:
        print(f"Could not find a solution for modulus q={q}")
        return None

    # Any non-zero column in the final basis is a solution
    final_solution = current_basis[:, 0]
    return final_solution

# Main execution
if __name__ == '__main__':
    # Problem parameters
    k_param = 3
    q_param = 2**k_param
    n_param = 4
    # Set m large enough to satisfy m > n*k, e.g. m = n*k + 1
    m_param = n_param * k_param + 1

    # Generate a random n x m matrix A with entries in Z_q
    np.random.seed(0)
    A_matrix = np.random.randint(0, q_param, size=(n_param, m_param))
    
    print("Finding a non-zero binary vector x such that Ax = 0 (mod q)")
    print(f"Parameters: n={n_param}, m={m_param}, q={q_param} (k={k_param})")
    
    # Solve the problem
    solution_x = solve_homogeneous_binary(A_matrix, q_param)

    if solution_x is not None:
        print("\nFound a solution vector x:")
        print(solution_x)

        # Verification
        result_vector = (A_matrix @ solution_x) % q_param
        print("\nVerification: A*x (mod q):")
        print(result_vector)

        is_correct = np.all(result_vector == 0)
        is_nonzero = np.any(solution_x != 0)
        is_binary = np.all((solution_x == 0) | (solution_x == 1))
        
        print(f"\nSolution is non-zero: {is_nonzero}")
        print(f"Solution is binary: {is_binary}")
        print(f"Solution is correct (Ax=0 mod q): {is_correct}")

        if is_correct and is_nonzero and is_binary:
            print("\nThe solution is valid.")
    else:
        print("\nNo solution was found by the algorithm.")

<<<E>>>