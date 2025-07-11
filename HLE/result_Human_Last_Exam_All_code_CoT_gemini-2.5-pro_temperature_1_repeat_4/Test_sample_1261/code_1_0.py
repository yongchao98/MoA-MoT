import numpy as np

def z2_null_space(A):
    """
    Computes the basis for the null space of a matrix A over Z_2.
    
    Args:
        A (np.ndarray): The input matrix with entries in {0, 1}.
        
    Returns:
        np.ndarray: A matrix whose columns form a basis for the null space of A.
                    Returns None if the null space is trivial.
    """
    mat = np.copy(A).astype(int)
    n_rows, n_cols = mat.shape
    pivot_row = 0
    pivot_cols = []

    for j in range(n_cols):
        if pivot_row < n_rows:
            pivot = np.where(mat[pivot_row:, j] == 1)[0]
            if pivot.size > 0:
                pivot_idx = pivot_row + pivot[0]
                mat[[pivot_row, pivot_idx]] = mat[[pivot_idx, pivot_row]]  # Swap rows
                
                for i in range(n_rows):
                    if i != pivot_row and mat[i, j] == 1:
                        mat[i, :] = (mat[i, :] + mat[pivot_row, :]) % 2
                
                pivot_cols.append(j)
                pivot_row += 1

    free_cols = [j for j in range(n_cols) if j not in pivot_cols]
    
    if not free_cols:
        return None

    dim_null_space = len(free_cols)
    null_basis = np.zeros((n_cols, dim_null_space), dtype=int)
    
    # Create basis vectors from free variables
    for i, free_col_idx in enumerate(free_cols):
        basis_vec = np.zeros(n_cols, dtype=int)
        basis_vec[free_col_idx] = 1
        for p_i, p_j in enumerate(pivot_cols):
            if mat[p_i, free_col_idx] == 1:
                basis_vec[p_j] = 1
        null_basis[:, i] = basis_vec

    return null_basis

def solve_homogeneous_mod_power_of_2(A, q):
    """
    Finds a non-zero vector x in {0,1}^m such that Ax = 0 (mod q),
    where q is a power of 2.

    Args:
        A (np.ndarray): The n x m matrix.
        q (int): The modulus, must be a power of 2.

    Returns:
        np.ndarray: A non-zero solution vector x in {0,1}^m.
    """
    if q == 0 or (q & (q - 1)) != 0:
        raise ValueError("q must be a power of 2.")
    
    k = int(np.log2(q))
    n, m = A.shape
    
    # Step 1: Find basis for solutions mod 2
    A_mod_2 = A % 2
    # Basis for L_1, columns are basis vectors
    basis = z2_null_space(A_mod_2)
    
    if basis is None:
        print("No non-trivial solution exists.")
        return None

    # Iteratively lift the solution from mod 2^j to mod 2^(j+1)
    for j in range(1, k):
        # For each basis vector, calculate Ax / 2^j
        d = basis.shape[1] # Dimension of the current solution space
        V = np.zeros((n, d), dtype=int)
        
        # This part can be slow if d is large, but shows the logic
        for i in range(d):
            b = basis[:, i]
            Ax_b = A @ b
            # Integer division, then mod 2
            v_i = (Ax_b // (2**j)) % 2
            V[:, i] = v_i

        # Find kernel of the transformation defined by V
        # i.e., find C such that V*C = 0 mod 2
        C = z2_null_space(V)

        if C is None:
            # This should not happen given the problem constraints (m > kn)
            print(f"Lifting failed at step j={j+1}. No further non-trivial solutions found.")
            # Return the best solution found so far (solves up to mod 2^j)
            return basis[:, 0]

        # Update the basis for the new solution space L_{j+1}
        # Basis for L_{j+1} is B_j * C
        basis = (basis @ C) % 2

    # The final basis is for the solution space mod q=2^k.
    # Any non-zero column of the final basis matrix is a valid solution.
    return basis[:, 0]

if __name__ == '__main__':
    # Set parameters for the problem
    # n=2, k=3 -> q=8. We need m > nk = 6. Let's choose m=7.
    n = 2
    k = 3
    q = 2**k
    m = 7
    
    # To ensure reproducibility, we use a fixed seed
    np.random.seed(42)
    # Generate a random matrix A with entries in Z_q
    A = np.random.randint(0, q, size=(n, m))
    
    print("Let q = 2^k.")
    print(f"Given parameters: n={n}, k={k}, q={q}, m={m}")
    print(f"The condition m > nk is {m} > {n*k}, which is {m > n*k}.")
    print("\nMatrix A (mod q):")
    print(A)

    # Find the solution
    solution_x = solve_homogeneous_mod_power_of_2(A, q)
    
    print("\nFound a non-zero solution x in {0,1}^m:")
    print(solution_x)
    
    # Verification
    if solution_x is not None:
        Ax = A @ solution_x
        Ax_mod_q = Ax % q
        
        print("\nVerification:")
        print(f"A * x = {Ax}")
        print(f"A * x (mod {q}) = {Ax_mod_q}")
        
        if np.all(Ax_mod_q == 0):
            print("The solution is correct. Ax is the zero vector mod q.")
        else:
            print("The solution is incorrect.")
        
        print("\nFinal Equation:")
        equation_str = []
        for i in range(m):
            if solution_x[i] == 1:
                equation_str.append(f"{solution_x[i]} * {A[:, i].tolist()}")
        print(" + ".join(equation_str) + f" = {Ax.tolist()} = 0 (mod {q})")

<<<E>>>