import numpy as np

def null_space_F2(A):
    """
    Finds a basis for the null space of a matrix A over F_2 using
    Reduced Row Echelon Form.
    Returns a matrix whose columns are the basis vectors of the null space.
    """
    rows, cols = A.shape
    M = np.copy(A).astype(int) % 2
    pivot_cols = []
    
    # Compute Reduced Row Echelon Form (RREF) over F2
    pivot_row = 0
    for j in range(cols):
        if pivot_row >= rows:
            break
        
        # Find a pivot (a '1') in the current column
        pivot = np.where(M[pivot_row:, j] == 1)[0]
        if pivot.size > 0:
            # Move pivot row to current position
            p_idx = pivot[0] + pivot_row
            M[[pivot_row, p_idx]] = M[[p_idx, pivot_row]]
            
            # Eliminate other '1's in the pivot column
            for i in range(rows):
                if i != pivot_row and M[i, j] == 1:
                    M[i] = (M[i] + M[pivot_row]) % 2
            
            pivot_cols.append(j)
            pivot_row += 1
            
    # RREF is computed. Now find the null space basis.
    rank = len(pivot_cols)
    free_cols = [j for j in range(cols) if j not in pivot_cols]
    
    basis = []
    for free_col in free_cols:
        # Create a basis vector for each free variable
        vec = np.zeros(cols, dtype=int)
        vec[free_col] = 1
        
        # Express pivot variables in terms of the free variable
        for i in range(rank):
            pivot_col = pivot_cols[i]
            if M[i, free_col] == 1:
                vec[pivot_col] = 1
        basis.append(vec)
        
    if not basis:
        return np.array([]).reshape((cols, 0))
    
    return np.vstack(basis).T

def solve_ax_zero_mod_power_of_2(A, k):
    """
    Finds a non-zero vector x in {0,1}^m such that Ax = 0 (mod 2^k).
    This implements the deterministic lifting algorithm.
    A is an n x m matrix.
    k is the exponent, so q = 2^k.
    """
    n, m = A.shape
    q = 2**k

    # Step 1: Solve Ax = 0 (mod 2)
    # Find a basis for the null space of A mod 2. The columns of X form the basis.
    A_mod_2 = A % 2
    X = null_space_F2(A_mod_2)
    
    if X.shape[1] == 0:
        print("Error: No non-trivial solution found mod 2.")
        print("This should not happen if m > n.")
        return None

    # Iteratively lift the solution from mod 2^i to mod 2^(i+1)
    for i in range(1, k):
        # The columns of the current matrix X are vectors that satisfy Ax = 0 (mod 2^i)
        
        # To lift to mod 2^(i+1), we need to find a linear combination of these vectors,
        # Xz, such that A(Xz) = 0 (mod 2^(i+1)).
        # Let v_j = Ax_j. We know v_j = 0 (mod 2^i), so v_j = 2^i * b_j.
        # We need sum(z_j * v_j) = 0 (mod 2^(i+1)), which means sum(z_j * b_j) = 0 (mod 2).
        # This is a new linear system Bz=0 (mod 2).
        
        # Form the matrix B for the new system.
        current_power_of_2 = 2**i
        V = (A @ X)
        B = (V // current_power_of_2) % 2

        # Find a basis Z for the null space of B.
        Z = null_space_F2(B)
        
        if Z.shape[1] == 0:
            # This means no non-trivial combination z exists. This implies one of the
            # original vectors in X was already a solution for the next level.
            # This case is unlikely given the problem constraints but we can handle it.
            print(f"Lifting found an existing solution at level {i+1}.")
            for j in range(X.shape[1]):
                if np.all(((A @ X[:, j]) % (2**(i+1))) == 0):
                    return X[:, j] # Return the found solution
            print("Error: Lifting failed unexpectedly.")
            return None

        # Update the basis for the next level. The new solutions are linear
        # combinations of the old solutions, with coefficients from Z.
        X = (X @ Z) % 2

    # After k-1 loops, any column in the final matrix X is a solution mod 2^k.
    # We just need one non-zero solution.
    solution_vector = X[:, 0]
    return solution_vector

def main():
    """Main function to demonstrate the solution."""
    # Setup problem parameters satisfying the prompt's constraints.
    n = 3
    k = 4  # k > 1
    # We need m = Ω(n^k), so m = Ω(3^4) = Ω(81). Let's pick m=85.
    # We also need m ∈ poly(n). m=85 is constant w.r.t n, so it's poly(n).
    m = 85
    
    q = 2**k

    # Generate a random matrix A, with a seed for reproducibility.
    np.random.seed(42)
    A = np.random.randint(0, q, size=(n, m))

    print("--- Problem Statement ---")
    print(f"Find a non-zero vector x in {{0,1}}^{m} such that Ax = 0 (mod q)")
    print(f"Parameters: n={n}, m={m}, k={k}, q={q}")
    # print(f"Matrix A (shape {n}x{m}):\n{A}\n") # Too large to print

    # Find the solution using the deterministic algorithm
    x = solve_ax_zero_mod_power_of_2(A, k)

    print("\n--- Algorithm Output ---")
    if x is not None:
        print(f"Found a non-zero binary solution vector x (first 20 elements shown):")
        print(f"{x[:20]}...\n")
        
        # Verification
        print("--- Verification ---")
        # The final equation is Ax = 0 (mod q). We print the vectors involved.
        # Printing A is too large. We will print the result of A @ x.
        print(f"Solution vector x (size {x.size}):\n{x}\n")

        result_vector = (A @ x) % q
        
        print(f"Result of (A @ x) mod {q}:")
        print(result_vector)

        if np.all(result_vector == 0):
            print("\nVerification successful: Ax is the zero vector mod q.")
        else:
            print("\nVerification failed: Ax is NOT the zero vector mod q.")
        
        if np.all(np.isin(x, [0, 1])) and np.any(x != 0):
            print("The vector x is indeed a non-zero binary vector.")
        else:
            print("The vector x is NOT a non-zero binary vector.")

    else:
        print("The algorithm did not find a solution.")

if __name__ == "__main__":
    main()