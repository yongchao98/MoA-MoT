import numpy as np
import sympy

def find_f2_null_space(M):
    """
    Finds the basis for the null space of a binary matrix M over F_2.
    Uses sympy's rref for robustness.
    Returns a matrix whose columns form the basis of the null space.
    """
    if M.shape[1] == 0:
        return np.empty((M.shape[1], 0), dtype=np.int64)
    
    # Use sympy to calculate the reduced row echelon form
    M_sympy = sympy.Matrix(M)
    # The mod=2 option ensures arithmetic is in F_2, but rref over integers
    # followed by mod 2 also works since we only have 0s and 1s.
    rref_M_sympy, pivot_cols = M_sympy.rref()
    rref_M = np.array(rref_M_sympy).astype(np.int64)
    
    n_vars = M.shape[1]
    free_cols = [i for i in range(n_vars) if i not in pivot_cols]
    
    basis = []
    for col_idx in free_cols:
        # Create a basis vector for each free variable
        vec = np.zeros(n_vars, dtype=np.int64)
        vec[col_idx] = 1
        # Back-substitute to find values for pivot variables
        for i, piv_col in enumerate(pivot_cols):
            if rref_M[i, col_idx] == 1:
                vec[piv_col] = 1  # In F_2, -1 is 1
        basis.append(vec)
        
    if not basis:
        return np.empty((n_vars, 0), dtype=np.int64)
        
    return np.array(basis, dtype=np.int64).T

def solve_for_binary_vector(A, q, k):
    """
    Implements the lifting algorithm to find a non-zero x in {0,1}^m
    such that Ax = 0 (mod q).
    """
    print("Starting lifting algorithm...")
    n, m = A.shape
    A_current = np.copy(A)
    bases = []

    # Lifting loop: from mod 2 up to mod 2^k = q
    for i in range(1, k + 1):
        print(f"\nStep {i}/{k}: Solving for mod {2**i}")
        print(f"Matrix A_{i-1} has shape {A_current.shape}")
        
        # System to solve is A_{i-1} * y = 0 (mod 2)
        M = A_current % 2
        
        # Find basis for the null space over F_2
        basis_matrix = find_f2_null_space(M)
        
        if basis_matrix.shape[1] == 0:
            print(f"Error: Null space at step {i} is trivial. Cannot proceed.")
            return None
        
        print(f"Found a null space of dimension {basis_matrix.shape[1]}.")
        bases.append(basis_matrix)
        
        # Update for the next step of lifting
        # A_i = (A_{i-1} * B_i) / 2
        if i < k:
            A_next = np.matmul(A_current, basis_matrix, dtype=np.int64)
            # All elements must be even by construction
            if np.any(A_next % 2 != 0):
                 print(f"Error: Matrix is not divisible by 2 at step {i}. Something went wrong.")
                 return None
            A_current = A_next // 2
    
    print("\nLifting complete. Reconstructing the solution vector x.")
    
    # Start with a non-zero vector from the final null space
    # (the first basis vector)
    x = bases[-1][:, 0]
    
    # Back-substitute through the bases
    for i in range(k - 2, -1, -1):
        x = np.matmul(bases[i], x, dtype=np.int64)
        
    return x

# --- Main execution ---
if __name__ == '__main__':
    # Setup parameters satisfying the problem constraints
    # k > 1
    # m = Omega(n^k) and m in poly(n)
    k = 2
    n = 4
    # Let m = n^k + 4 = 4^2 + 4 = 20
    # This satisfies m > n*k (20 > 8)
    m = n**k + 4
    q = 2**k

    print(f"Problem setup:")
    print(f"q = 2^k = {q} (k={k})")
    print(f"Matrix A is {n}x{m}")
    print(f"Constraint check: m={m}, n^k={n**k}, nk={n*k}. m > nk is satisfied.")
    
    # Generate a random matrix A
    A = np.random.randint(0, q, size=(n, m), dtype=np.int64)
    
    # Solve the problem
    x = solve_for_binary_vector(A, q, k)

    # Print results
    if x is not None:
        print("\n--- Algorithm Finished ---")
        print("Found a non-zero binary vector x.")
        
        # Verify the solution
        is_binary = np.all((x == 0) | (x == 1))
        is_nonzero = np.any(x != 0)
        result_vector = np.matmul(A, x, dtype=np.int64) % q
        is_solution = np.all(result_vector == 0)

        print("\nVerification:")
        print(f"Is x binary? {is_binary}")
        print(f"Is x non-zero? {is_nonzero}")
        print(f"Is Ax = 0 (mod q)? {is_solution}")

        print("\n--- Final Equation ---")
        # To avoid printing a huge matrix, we show a snippet
        # and print the vectors.
        print(f"q = {q}")
        print(f"A (first 5x5 snippet):\n{A[:5, :5]}\n...")
        print(f"x =\n{x}")
        print(f"Result (Ax mod q) =\n{result_vector}")
        
        if not is_solution:
            print("\nError: The found vector is not a correct solution.")
