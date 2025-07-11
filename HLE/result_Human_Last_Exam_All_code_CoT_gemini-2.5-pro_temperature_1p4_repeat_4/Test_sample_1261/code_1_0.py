import numpy as np

def find_null_space_basis_mod2(A):
    """
    Finds a basis for the null space of a matrix A over F_2 using Gaussian elimination.

    Args:
        A (np.ndarray): The input matrix with integer entries.

    Returns:
        list[np.ndarray]: A list of basis vectors for the null space.
                          Returns an empty list if only the trivial solution exists.
    """
    # Work with a copy of the matrix modulo 2
    mat = A.copy() % 2
    n, m = mat.shape
    
    pivot_row = 0
    pivot_cols = []
    
    # Forward elimination to get row echelon form
    for col in range(m):
        if pivot_row < n:
            pivot = np.where(mat[pivot_row:, col] == 1)[0]
            if len(pivot) > 0:
                # Swap rows to bring pivot to the current pivot_row
                pivot_index_in_slice = pivot[0]
                actual_pivot_row = pivot_row + pivot_index_in_slice
                mat[[pivot_row, actual_pivot_row]] = mat[[actual_pivot_row, pivot_row]]
                
                # Eliminate other 1s in the same column
                for i in range(n):
                    if i != pivot_row and mat[i, col] == 1:
                        mat[i] = (mat[i] + mat[pivot_row]) % 2
                
                pivot_cols.append(col)
                pivot_row += 1

    # Identify free columns
    free_cols = [j for j in range(m) if j not in pivot_cols]
    
    # Back substitution to find basis vectors
    basis = []
    for free_col in free_cols:
        # Create a vector for this free variable
        x = np.zeros(m, dtype=int)
        x[free_col] = 1
        
        # Solve for pivot variables
        for i in range(len(pivot_cols) - 1, -1, -1):
            pivot_c = pivot_cols[i]
            pivot_r = i # In reduced form, pivot_row i has pivot in col pivot_cols[i]
            
            # The value of the pivot variable is the dot product of its row
            # (excluding the pivot element itself) and the solution vector x.
            val = np.dot(mat[pivot_r, :], x) % 2
            x[pivot_c] = val
            
        basis.append(x)
        
    return basis

def solve_short_integer_solution(A, q):
    """
    Finds a non-zero binary vector x such that Ax = 0 (mod q),
    where q is a power of 2.

    Args:
        A (np.ndarray): An n x m matrix with entries in Z_q.
        q (int): The modulus, must be a power of 2, q = 2^k.

    Returns:
        np.ndarray: A non-zero binary vector x of length m, or None if failed.
    """
    if not (q > 0 and (q & (q - 1)) == 0):
        raise ValueError("q must be a power of 2.")
    if q == 1:
        # Trivial case, any non-zero x works. Let's return one.
        x = np.zeros(A.shape[1], dtype=int)
        x[0] = 1
        return x

    k = int(np.log2(q))
    n, m = A.shape

    # `current_cols` stores vectors that are sums of original columns of A
    # `tracker_vectors` stores how each current_col is formed from the original columns
    current_cols = [A[:, j] for j in range(m)]
    tracker_vectors = [np.eye(m, dtype=int)[j] for j in range(m)]

    for i in range(1, k + 1):
        # At the start of this loop, all vectors in `current_cols` are 0 (mod 2**(i-1))
        
        # Create the matrix for the mod 2 system
        # The division by 2**(i-1) is exact integer division
        mod2_matrix = np.array([col // (2**(i-1)) for col in current_cols]).T
        
        # Find basis for the null space mod 2
        null_space_basis = find_null_space_basis_mod2(mod2_matrix)

        if not null_space_basis:
            print(f"Error: Could not find a non-trivial solution at lifting step {i}.")
            print(f"This should not happen if m > n*k.")
            return None

        # Prepare for the next iteration
        next_cols = []
        next_trackers = []
        
        for basis_vec in null_space_basis:
            new_col = np.zeros(n, dtype=int)
            new_tracker = np.zeros(m, dtype=int)

            for j in range(len(basis_vec)):
                if basis_vec[j] == 1:
                    new_col = new_col + current_cols[j]
                    new_tracker = (new_tracker + tracker_vectors[j]) % 2
            
            next_cols.append(new_col)
            next_trackers.append(new_tracker)
            
        current_cols = next_cols
        tracker_vectors = next_trackers
    
    # After k loops, any vector in current_cols is 0 (mod 2**k)
    # The corresponding tracker vector is the solution x
    # We just need one non-zero solution.
    for x in tracker_vectors:
        if np.any(x):
            return x

    return None

if __name__ == '__main__':
    # Setup parameters satisfying the problem constraints
    # Let n=2, k=3. Then q = 8.
    # The constraint m = Omega(n^k) = Omega(2^3=8) means m must be large.
    # The iterative algorithm requires m > n*k = 2*3 = 6 to guarantee a solution.
    # We choose m that satisfies this simpler condition for a demonstration.
    n = 2
    k = 3
    q = 2**k
    m = n * k + 2  # m=8, which is > n*k=6

    print(f"Parameters: n={n}, m={m}, k={k}, q={q}")
    print(f"Constraint check: m > n*k is {m > n*k}")
    
    # Generate a random n x m matrix A with entries in Z_q
    A = np.random.randint(0, q, size=(n, m))
    print("\nGenerated random matrix A:")
    print(A)

    # Solve for a non-zero binary x such that Ax = 0 (mod q)
    x = solve_short_integer_solution(A, q)

    if x is not None:
        print("\nFound a non-zero binary solution x:")
        print(x)

        # Verification
        result_vector = A @ x
        print("\nVerification:")
        print("A @ x = ")
        print(result_vector)
        print(f"A @ x (mod {q}) = ")
        print(result_vector % q)

        # Final check that the result is indeed the zero vector mod q
        if np.all((A @ x) % q == 0):
            print("\nSuccess: The solution is correct.")
        else:
            print("\nFailure: The solution is incorrect.")
    else:
        print("\nCould not find a solution.")
