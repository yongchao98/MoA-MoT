import numpy as np

def find_f2_nullspace_basis(M):
    """
    Finds a basis for the null space of a matrix M over F_2 using Gaussian elimination.
    
    Args:
        M: A numpy array with entries in {0, 1}.
        
    Returns:
        A list of numpy arrays, where each array is a basis vector of the null space.
    """
    n_rows, n_cols = M.shape
    M_aug = np.copy(M)
    
    pivot_cols = []
    pivot_row = 0
    for j in range(n_cols):
        if pivot_row < n_rows and M_aug[pivot_row, j] == 1:
            # Pivot is already in place
            pass
        elif pivot_row < n_rows:
            # Search for a row to swap
            found = False
            for i in range(pivot_row + 1, n_rows):
                if M_aug[i, j] == 1:
                    M_aug[[pivot_row, i]] = M_aug[[i, pivot_row]] # Swap rows
                    found = True
                    break
            if not found:
                continue # No pivot in this column

        # Eliminate other 1s in the current column
        for i in range(n_rows):
            if i != pivot_row and M_aug[i, j] == 1:
                M_aug[i, :] = (M_aug[i, :] + M_aug[pivot_row, :]) % 2

        pivot_cols.append(j)
        pivot_row += 1
        if pivot_row >= n_rows:
            break
            
    # Back substitution to find basis vectors
    free_cols = [j for j in range(n_cols) if j not in pivot_cols]
    basis = []
    
    for free_col in free_cols:
        b = np.zeros(n_cols, dtype=int)
        b[free_col] = 1
        for i in range(len(pivot_cols) - 1, -1, -1):
            pivot_c = pivot_cols[i]
            row_sum = np.dot(M_aug[i, pivot_c+1:], b[pivot_c+1:])
            b[pivot_c] = row_sum % 2
        basis.append(b)
        
    return basis

def solve_for_binary_null_vector(A, q, k):
    """
    Finds a non-zero vector x in {0,1}^m such that Ax = 0 (mod q),
    where q = 2^k.
    """
    n, m = A.shape
    
    # Check if a solution is guaranteed by the lifting method
    if m <= n * k:
        print(f"Warning: m={m} is not guaranteed to be large enough "
              f"compared to n*k={n*k}. Algorithm may fail.")

    # S_i contains vectors x s.t. Ax = 0 (mod 2^i)
    # Start with S_0 = {e_1, ..., e_m}, the standard basis vectors for R^m.
    # Ax = 0 (mod 2^0=1) is trivially true for all x.
    candidate_vectors = [v for v in np.eye(m, dtype=int)]
    
    divisor = 1
    for i in range(1, k + 1):
        # We have vectors x in candidate_vectors s.t. Ax = 0 (mod divisor).
        # So Ax = divisor * v for integer vectors v.
        # We want to find a combination sum(c_j * x_j) s.t.
        # A * sum(c_j * x_j) = 0 (mod divisor*2).
        # This means sum(c_j * (Ax_j / divisor)) = 0 (mod 2).
        
        M_cols = []
        for x_vec in candidate_vectors:
            v = np.dot(A, x_vec)
            # This division is guaranteed to be exact.
            v_prime = v // divisor
            M_cols.append(v_prime)
        
        if not M_cols:
            print(f"Algorithm failed at step {i}: no candidate vectors left.")
            return None

        M = np.array(M_cols).T
        M_f2 = M % 2
        
        null_basis = find_f2_nullspace_basis(M_f2)

        if not null_basis:
            print(f"Algorithm failed at step {i}: no non-trivial solution found for mod 2 system.")
            return None

        # Build the next set of candidate vectors
        new_candidates = []
        for b in null_basis:
            # b is a basis vector from the nullspace of M_f2.
            # It tells us which of the current candidates to sum up.
            new_x = np.zeros(m, dtype=int)
            for j in range(len(b)):
                if b[j] == 1:
                    new_x = (new_x + np.array(candidate_vectors[j]))
            # The coefficients are summed over integers, then taken mod 2 for the final {0,1} vector
            new_x = new_x % 2
            
            # Ensure the solution is non-zero
            if np.any(new_x):
                new_candidates.append(new_x)
        
        candidate_vectors = new_candidates
        divisor *= 2

    if not candidate_vectors:
        print("Algorithm terminated without finding a solution.")
        return None
    
    # Any vector in candidate_vectors is a solution
    return candidate_vectors[0]

def main():
    # Setup parameters satisfying the problem constraints
    # Let n=2, k=3. Then q = 2^3 = 8.
    # We need m = Omega(n^k) = Omega(2^3=8) and m > n*k = 2*3=6.
    # Let's choose m=8.
    n = 2
    k = 3
    m = 8
    q = 2**k

    print(f"Parameters: n={n}, m={m}, k={k}, q={q}")
    print("-" * 30)

    # Generate a random matrix A from Z_q^(n x m)
    # For reproducibility, we use a fixed seed
    np.random.seed(0)
    A = np.random.randint(0, q, size=(n, m))
    
    print("Input Matrix A (mod q):")
    print(A)
    print("-" * 30)
    
    # Find the solution vector x
    x = solve_for_binary_null_vector(A, q, k)

    if x is not None:
        print("Found a non-zero solution vector x in {0,1}^m:")
        print(x)
        print("-" * 30)

        # Verification
        result_vector = np.dot(A, x)
        
        print("Verification: Ax (mod q)")
        print(f"A =\n{A}\n")
        print(f"x =\n{x.reshape(-1, 1)}\n")
        print(f"Ax =\n{result_vector.reshape(-1, 1)}\n")
        print(f"Ax mod {q} =\n{(result_vector % q).reshape(-1, 1)}")

        # Final check
        if np.all((result_vector % q) == 0):
            print("\nSuccess: Ax is indeed 0 (mod q).")
        else:
            print("\nFailure: Ax is not 0 (mod q).")
    else:
        print("Could not find a solution.")

if __name__ == "__main__":
    main()