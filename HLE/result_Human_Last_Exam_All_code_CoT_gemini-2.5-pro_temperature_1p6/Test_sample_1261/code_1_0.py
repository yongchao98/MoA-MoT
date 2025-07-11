import numpy as np

def gaussian_elimination_mod2(matrix):
    """
    Finds a basis for the null space of a matrix over F_2 using Gaussian elimination.
    Input: A numpy array with integer entries.
    Output: A list of numpy arrays, representing the basis vectors for the null space.
            Returns an empty list if only the trivial solution exists.
    """
    if matrix.size == 0:
        return []
    A = matrix.copy() % 2
    n, m = A.shape
    pivot_row, pivot_col = 0, 0
    pivot_cols = []
    
    # Create a copy for row reduction
    reduced_A = A.copy()

    # Forward elimination to get row echelon form
    while pivot_row < n and pivot_col < m:
        i = pivot_row
        while i < n and reduced_A[i, pivot_col] == 0:
            i += 1
        
        if i < n:  # Pivot found
            reduced_A[[pivot_row, i]] = reduced_A[[i, pivot_row]]  # Swap rows
            # Eliminate other 1s in the pivot column for all rows
            for j in range(n):
                if j != pivot_row and reduced_A[j, pivot_col] == 1:
                    reduced_A[j, :] = (reduced_A[j, :] + reduced_A[pivot_row, :]) % 2
            pivot_cols.append(pivot_col)
            pivot_row += 1
        pivot_col += 1
        
    # Find null space basis from the reduced form
    null_space_basis = []
    free_cols = [j for j in range(m) if j not in pivot_cols]
    
    for free_col in free_cols:
        x = np.zeros(m, dtype=int)
        x[free_col] = 1
        for i, p_col in enumerate(pivot_cols):
            x[p_col] = reduced_A[i, free_col]
        null_space_basis.append(x)
        
    return null_space_basis

def solve_binary_sis_with_lifting(A, n, m, k):
    """
    Solves Ax = 0 (mod 2^k) for a binary vector x using the lifting algorithm.
    """
    q = 2**k
    
    print(f"Solving Ax = 0 (mod {q}) for A in Z_{q}^{{{n}x{m}}}")
    print(f"Goal: find a non-zero x in {{0,1}}^{m}")
    print("=" * 40)
    print("Algorithm: Deterministic lifting from mod 2 to mod 2^k.\n")

    # The condition m > k*n ensures the algorithm's success.
    # The problem statement's constraints guarantee this for large n.
    if m <= n * k:
        print(f"Warning: Condition m > kn not met (m={m}, kn={k*n}). Success is not guaranteed.")

    # --- Step 1: Base case, solve Ax = 0 (mod 2) ---
    print("--- Step 1: Solving Ax = 0 (mod 2) ---")
    current_basis = gaussian_elimination_mod2(A)
    
    if not current_basis:
        print("Failed: No non-trivial solution for Ax = 0 (mod 2). This shouldn't happen if m > n.")
        return None
    
    print(f"Found a basis of dimension {len(current_basis)} for the solution space mod 2.\n")
    
    # --- Step 2: Iteratively lift the solution from mod 2^i to mod 2^(i+1) ---
    for i in range(1, k):
        current_mod = 2**i
        next_mod = 2**(i+1)
        print(f"--- Lifting Step {i}: from mod {current_mod} to mod {next_mod} ---")
        
        d_i = len(current_basis)
        print(f"Current solution space dimension: {d_i}")

        if d_i == 0:
             print("Lifting failed: solution space is empty.")
             return None

        # Form the lifting matrix Y
        Y_cols = []
        for b_vec in current_basis:
            # By induction, A @ b_vec is divisible by current_mod
            y_vec = (A @ b_vec) // current_mod
            Y_cols.append(y_vec)
        
        # Y is an n x d_i matrix
        Y = np.array(Y_cols).T
        
        # Find a non-trivial vector `c` in the null space of Y (mod 2)
        print(f"Solving Yc = 0 (mod 2) for lifting matrix Y of size {Y.shape}")
        c_basis = gaussian_elimination_mod2(Y)
        
        if not c_basis:
            print("Lifting failed: no non-trivial solution for Yc = 0 (mod 2).")
            print("This can happen if the dimension of the solution space becomes too small.")
            return None
        
        # Construct the new basis for the solution space mod `next_mod`
        new_basis = []
        current_basis_mat = np.array(current_basis).T
        for c_vec in c_basis:
            # The new basis vectors are linear combinations of the old ones.
            # new_b = sum(c_j * old_b_j) mod 2
            new_b = (current_basis_mat @ c_vec) % 2
            new_basis.append(new_b)
        
        current_basis = new_basis
        print(f"Successfully lifted. New solution space dimension: {len(current_basis)}\n")

    # After the loop, any non-zero vector in current_basis is a solution to Ax=0 (mod 2^k)
    if not current_basis:
        print("Algorithm finished but found no non-trivial solution.")
        return None

    final_solution = current_basis[0]

    # --- Step 3: Output and verify the final solution ---
    print("=" * 40)
    print("Algorithm finished. Found a solution.")
    print("Final non-zero solution vector x:")
    print(final_solution)
    
    # Verify the solution
    print("\n--- Verification ---")
    Ax = A @ final_solution
    Ax_mod_q = Ax % q
    
    print(f"A * x = {Ax}")
    print(f"A * x (mod {q}) = {Ax_mod_q}")
    
    if np.all(Ax_mod_q == 0):
        print("Verification SUCCESS: Ax = 0 (mod q).")
    else:
        print("Verification FAILED.")

    print("\n--- Final Equation Rendered ---")
    # Output the equation in the format: c1*v1 + c2*v2 + ... = result (mod q)
    terms = []
    for i in range(m):
        # Python's default string conversion for numpy arrays is clear enough
        terms.append(f"{final_solution[i]} * {A[:, i]}")
    equation = " + ".join(terms)
    print(f"{equation} = {Ax} \nWhich is {Ax_mod_q} (mod {q})")

if __name__ == '__main__':
    # Set up parameters according to the problem statement.
    n_param = 3  # n
    k_param = 4  # k > 1
    # m = Ω(n^k). m must be > n*k for the algorithm to be guaranteed success.
    m_param = n_param * k_param + 2 # Let m = 3*4 + 2 = 14
    q_param = 2**k_param

    # Generate a random matrix A for demonstration
    np.random.seed(42)
    A_matrix = np.random.randint(0, q_param, size=(n_param, m_param))
    
    solve_binary_sis_with_lifting(A_matrix, n_param, m_param, k_param)