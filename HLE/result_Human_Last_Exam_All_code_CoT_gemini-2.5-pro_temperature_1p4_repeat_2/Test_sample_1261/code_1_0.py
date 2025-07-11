import numpy as np

def find_null_space_gf2(A):
    """
    Finds a basis for the null space of matrix A over GF(2) using Gaussian elimination.
    The matrix A is an n x m matrix. The function returns an m x d matrix, where d is
    the dimension of the null space, and its columns form the basis.
    """
    n, m = A.shape
    M = np.copy(A)
    
    pivot_row = 0
    pivot_cols = []
    # Convert to Reduced Row Echelon Form (RREF)
    for j in range(m):  # Iterate through columns
        if pivot_row >= n:
            break
        pivot_idx = -1
        for i in range(pivot_row, n):
            if M[i, j] == 1:
                pivot_idx = i
                break
        
        if pivot_idx != -1:
            # Swap rows to bring pivot to the top of the remaining submatrix
            M[[pivot_row, pivot_idx]] = M[[pivot_idx, pivot_row]]
            
            # Eliminate other 1s in the pivot column
            for i in range(n):
                if i != pivot_row and M[i, j] == 1:
                    M[i] = (M[i] + M[pivot_row]) % 2
            
            pivot_cols.append(j)
            pivot_row += 1

    # RREF is now computed in M. Identify free columns.
    free_cols = [j for j in range(m) if j not in pivot_cols]
    
    # Construct null space basis vectors
    basis = []
    for free_col in free_cols:
        v = np.zeros(m, dtype=int)
        v[free_col] = 1
        for i, pivot_col in enumerate(pivot_cols):
            if M[i, free_col] == 1:
                v[pivot_col] = 1
        basis.append(v)
    
    if not basis:
        return np.zeros((m, 0), dtype=int) # Empty basis for trivial null space
        
    return np.array(basis).T

def solve_binary_null_space(n, m, k):
    """
    Implements the lifting algorithm to find a non-zero x in {0,1}^m
    such that A*x = 0 (mod 2^k).
    
    The algorithm requires m > n*k for a guaranteed non-zero solution.
    """
    if m <= n * k:
        print(f"Warning: The condition m > n*k is not met ({m} <= {n*k}).")
        print("The algorithm is not guaranteed to find a non-zero solution.")

    q = 2**k
    # Set a seed for reproducibility
    np.random.seed(42)
    # Generate a random n x m matrix A with entries in Z_q
    A = np.random.randint(0, q, size=(n, m))
    
    print(f"Solving for a non-zero binary vector x such that Ax = 0 (mod {q})")
    print(f"Parameters: n={n}, m={m}, k={k}")

    # The lifting process starts with the full space Z_2^m, so the initial basis
    # for the solution space (mod 2^0=1) is the identity matrix.
    # The columns of `basis_B` are the basis vectors for the solution space.
    basis_B = np.identity(m, dtype=np.int64)

    # Loop from i=1 to k to lift the solution from mod 2^(i-1) to mod 2^i
    for i in range(1, k + 1):
        # Current dimension of the solution space
        d = basis_B.shape[1]
        
        if d == 0:
            print(f"Solution space became trivial at step i={i}. No non-zero solution found.")
            return

        # Define the linear map L(x) = (Ax / 2^(i-1)) mod 2.
        # We compute the matrix M for this map where columns are L(b) for b in basis_B.
        # Using np.int64 for matrix products to avoid overflow for reasonable inputs.
        prod = A.astype(np.int64) @ basis_B
        
        M_i_num = prod // (2**(i - 1))
        M_i = M_i_num % 2
        
        # Find the null space of M_i. The basis vectors (columns of C) tell us
        # how to combine vectors in basis_B to get solutions for the next level.
        C = find_null_space_gf2(M_i)
        
        # Update the basis for the solution space
        basis_B = (basis_B @ C) % 2
        
    final_dim = basis_B.shape[1]
    
    if final_dim == 0:
        print("\nFinal solution space is trivial. No non-zero binary solution found.")
        return

    # Any non-zero column of the final basis_B is a valid solution.
    # We pick the first one.
    x = basis_B[:, 0]
    
    print("\n--- RESULT ---")
    print("Found a non-zero binary vector x.")
    
    # Print the equation and verification
    print(f"\nThe equation is A * x = 0 (mod q), with q = {q}")
    print("\nMatrix A:")
    print(A)
    print("\nSolution vector x:")
    print(x)
    
    result_vector = (A @ x) % q
    print("\nVerification: A * x (mod q) results in:")
    print(result_vector)
    
    if np.all(result_vector == 0):
      print("\nVerification successful.")
    else:
      print("\nVerification FAILED.")

    print("\nThis vector x corresponds to a subset sum of columns of A that equals 0.")
    indices = np.where(x == 1)[0]
    print(f"The sum of columns at indices {list(indices)} of A should be the zero vector.")
    
    # Show the sum for each row
    for i in range(n):
        sum_val = np.sum(A[i, indices])
        lhs_parts = [str(A[i, j]) for j in indices]
        lhs = " + ".join(lhs_parts) if lhs_parts else "0"
        print(f"  Row {i}: {lhs} = {sum_val}, and {sum_val} mod {q} = {sum_val % q}")


# Example usage with parameters satisfying m > n*k
# Let n=3, k=2 -> q=4. We need m > 3*2=6. Let's use m=7.
solve_binary_null_space(n=3, m=7, k=2)
