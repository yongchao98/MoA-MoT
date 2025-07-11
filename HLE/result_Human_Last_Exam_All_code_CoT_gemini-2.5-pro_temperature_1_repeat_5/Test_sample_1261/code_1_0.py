import numpy as np

def find_kernel_z2(M):
    """
    Finds a basis for the kernel (null space) of a matrix M over Z_2 using Gaussian elimination.
    
    Args:
        M (np.ndarray): The input matrix with entries in {0, 1}.
    
    Returns:
        list of np.ndarray: A list of basis vectors for the kernel of M.
    """
    n, m = M.shape
    A = np.copy(M)
    
    pivot_col = 0
    pivot_rows = []
    
    # Forward elimination to get row echelon form
    for i in range(n):
        if pivot_col >= m:
            break
        
        # Find a pivot in the current column
        pivot_row = i
        while pivot_row < n and A[pivot_row, pivot_col] == 0:
            pivot_row += 1
        
        if pivot_row < n:
            A[[i, pivot_row]] = A[[pivot_row, i]] # Swap rows
            # Eliminate other 1s in the pivot column
            for row in range(n):
                if row != i and A[row, pivot_col] == 1:
                    A[row, :] = (A[row, :] + A[i, :]) % 2
            pivot_rows.append(i)
            pivot_col += 1
        else:
            # No pivot in this column, move to the next
            pivot_col += 1

    # Identify pivot and free columns
    rank = len(pivot_rows)
    pivot_cols = []
    col_idx = 0
    for row_idx in pivot_rows:
        while col_idx < m and A[row_idx, col_idx] == 0:
            col_idx += 1
        if col_idx < m:
            pivot_cols.append(col_idx)
        col_idx += 1
        
    free_cols = [j for j in range(m) if j not in pivot_cols]
    
    # Construct kernel basis vectors from free variables
    kernel_basis = []
    for free_col in free_cols:
        sol = np.zeros(m, dtype=int)
        sol[free_col] = 1
        for i, p_col in enumerate(pivot_cols):
            sol[p_col] = A[pivot_rows[i], free_col]
        kernel_basis.append(sol)
        
    return kernel_basis

def solve_subset_sum_mod_2k(A, k):
    """
    Finds a non-zero x in {0,1}^m such that Ax = 0 (mod 2^k).
    
    Args:
        A (np.ndarray): The n x m input matrix.
        k (int): The exponent for the modulus q=2^k.
        
    Returns:
        np.ndarray or None: A solution vector x, or None if no non-zero solution is found.
    """
    n, m = A.shape
    
    # Initially, the solution space mod 2^0=1 is the entire space {0,1}^m.
    # A basis for this space is the set of standard basis vectors.
    basis = [v for v in np.eye(m, dtype=int)]

    # Iteratively lift the solution from mod 2^(j-1) to mod 2^j
    for j in range(1, k + 1):
        d = len(basis)
        if d == 0:
            return None # No non-zero solutions exist
            
        # For each basis vector y_i of the current solution space, calculate
        # z_i = (A @ y_i) / 2^(j-1).
        # We then form a new system M*c = 0 (mod 2) where M has columns z_i mod 2.
        M = np.zeros((n, d), dtype=int)
        power_of_2 = 2**(j - 1)
        for i, y_i in enumerate(basis):
            # By induction, A @ y_i is divisible by power_of_2
            z_i = (A @ y_i) // power_of_2
            M[:, i] = z_i % 2
        
        # Find the kernel of M. The coefficients c will form the new basis.
        kernel_of_M = find_kernel_z2(M)
        
        if not kernel_of_M:
            return None # Kernel is trivial
            
        # The new basis for solutions mod 2^j is formed by linear combinations
        # of the old basis vectors, using the kernel vectors as coefficients.
        new_basis = []
        for c in kernel_of_M:
            new_y = np.zeros(m, dtype=int)
            for i in range(d):
                if c[i] == 1:
                    new_y = (new_y + basis[i]) % 2
            # Only add non-zero vectors to the new basis
            if np.any(new_y):
                 new_basis.append(new_y)
        basis = new_basis
    
    # After k loops, 'basis' is a basis for the solution space of Ax=0 (mod 2^k).
    # Return the first non-zero basis vector.
    return basis[0] if basis else None

if __name__ == '__main__':
    # Set parameters for the problem
    # n=2, k=2 -> q=4. We need m = Omega(n^k) = Omega(4). Let's choose m=6.
    n = 2
    k = 2
    m = 6
    q = 2**k

    # Ensure parameters satisfy the condition for the algorithm to likely succeed
    if m <= n * k:
        print(f"Warning: m={m} may not be large enough compared to n*k={n*k}.")
        print("The algorithm requires m to be large enough for the solution space to remain non-trivial.")

    # Generate a random matrix A
    np.random.seed(42) # for reproducibility
    A = np.random.randint(0, q, size=(n, m))

    print(f"Solving Ax = 0 (mod {q}) for x in {{0,1}}^{m}")
    print("------------------------------------------")
    print(f"Parameters: n={n}, m={m}, k={k}")
    print("\nRandomly sampled matrix A:")
    print(A)

    # Solve the problem
    x = solve_subset_sum_mod_2k(A, k)

    if x is not None:
        print("\nFound a non-zero solution x:")
        print(x)
        
        # Verification
        res = (A @ x)
        res_mod_q = res % q
        
        print("\n--- Verification Equation ---")
        for i in range(n):
            line = " ".join([f"{A[i,j]}*{x[j]}" for j in range(m)])
            print(f"({line}) = {res[i]}, which is {res_mod_q[i]} (mod {q})")
        
        print("\nFinal result vector Ax (mod q):")
        print(res_mod_q)
        
        is_solution = np.all(res_mod_q == 0)
        print(f"\nIs Ax = 0 (mod {q})? {is_solution}")
    else:
        print("\nNo non-zero binary solution was found.")