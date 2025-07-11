import numpy as np

def solve_binary_hsis(A, q):
    """
    Solves Ax = 0 (mod q) for a non-zero binary vector x using modulus lifting.
    
    Args:
        A (numpy.ndarray): The n x m input matrix with integer entries.
        q (int): The modulus, must be a power of 2 (q = 2^k, k > 1).

    Returns:
        numpy.ndarray: A non-zero binary vector x of size m, or None if no solution is found.
    """
    
    # --- Setup and Input Validation ---
    n, m = A.shape
    if q <= 1 or (q & (q - 1)) != 0:
        print(f"Error: Modulus q must be a power of 2 greater than 1. Got q={q}.")
        return None
    k = int(np.log2(q))

    print("--- Problem Setup ---")
    print(f"Goal: Find non-zero x in {{0,1}}^{m} such that Ax = 0 (mod {q})")
    print(f"Matrix A is {n}x{m}, q = 2^{k} = {q}.")
    
    # The condition m > nk ensures the algorithm's success.
    if m <= n * k:
        print(f"\nWarning: The theoretical guarantee requires m > n*k.")
        print(f"Current values: m={m}, n*k={n*k}. The algorithm may fail.")

    def null_space_f2(M):
        """
        Computes a basis for the null space of a matrix M over F_2.
        Returns a matrix whose columns are the basis vectors.
        """
        mat = np.copy(M) % 2
        rows, cols = mat.shape
        pivot_cols = []
        pivot_row_idx = 0
        for j in range(cols): # Iterate through columns
            if pivot_row_idx >= rows:
                break
            pivot_row = pivot_row_idx
            # Find a row with a 1 in this column to be the pivot
            while pivot_row < rows and mat[pivot_row, j] == 0:
                pivot_row += 1

            if pivot_row < rows:  # Pivot found
                mat[[pivot_row_idx, pivot_row]] = mat[[pivot_row, pivot_row_idx]] # Swap rows
                # Eliminate other 1s in this column
                for i in range(rows):
                    if i != pivot_row_idx and mat[i, j] == 1:
                        mat[i, :] = (mat[i, :] + mat[pivot_row_idx, :]) % 2
                pivot_cols.append(j)
                pivot_row_idx += 1
        
        free_cols = [j for j in range(cols) if j not in pivot_cols]
        basis = []

        for free_j in free_cols:
            b = np.zeros(cols, dtype=int)
            b[free_j] = 1
            # Solve for pivot variables
            for i, pivot_j in enumerate(pivot_cols):
                if mat[i, free_j] == 1:
                    b[pivot_j] = 1
            basis.append(b)

        if not basis:
            return np.zeros((cols, 0), dtype=int)
            
        return np.array(basis).T

    # --- Lifting Algorithm ---
    print("\n--- Starting Modulus Lifting ---")
    
    # Step 1: Base case for mod 2
    print("Step 1: Finding basis for solutions mod 2.")
    A_mod_2 = A % 2
    # B is a matrix where columns form the basis for the solution space.
    B = null_space_f2(A_mod_2)
    
    if B.shape[1] == 0:
        print("Error: No non-trivial binary solution found for mod 2. Cannot proceed.")
        return None
    print(f"Found basis of size {B.shape[1]} for solutions mod 2.")

    # Step 2: Iteratively lift from mod 2^j to mod 2^(j+1)
    for j in range(1, k):
        mod_current = 2**j
        mod_next = 2**(j+1)
        print(f"\nStep {j+1}: Lifting from mod {mod_current} to mod {mod_next}.")

        # V = A * B (over integers)
        V = np.dot(A, B)
        
        # U = (V / mod_current) mod 2
        U = (V // mod_current) % 2
        
        # Find basis for nullspace of U. C's columns are coefficients for combining B's columns.
        C = null_space_f2(U)
        
        if C.shape[1] == 0:
            print(f"Error: Solution space collapsed at this step. No non-zero solution found for mod {mod_next}.")
            return None
            
        # Update the basis B for the next iteration
        # New basis vectors are combinations of old ones.
        B = np.dot(B, C) % 2

        print(f"Found new basis of size {B.shape[1]} for solutions mod {mod_next}.")
    
    # --- Final Solution and Verification ---
    print("\n--- Algorithm Finished ---")
    if B.shape[1] == 0:
        print("Failed to find a non-zero solution.")
        return None

    # Any column of the final matrix B is a valid solution.
    x = B[:, 0]
    
    print(f"\nFound a non-zero binary solution vector x:\nx = {x.tolist()}")

    res = np.dot(A, x)
    res_mod_q = res % q
    
    print("\n--- Verification ---")
    print(f"A * x = {res.tolist()}")
    print(f"A * x (mod {q}) = {res_mod_q.tolist()}")

    if np.all(res_mod_q == 0):
        print("Verification successful: Ax is indeed 0 (mod q).")
    else:
        print("Verification FAILED: Ax is not 0 (mod q).")

    print("\n--- Final Equation Details ---")
    for i in range(n):
        eq_str_parts = []
        for j_idx in range(m):
            # Only include terms where x_j is 1
            if x[j_idx] == 1:
                eq_str_parts.append(f"{A[i, j_idx]}*{x[j_idx]}")
        eq_str = " + ".join(eq_str_parts)
        print(f"Row {i+1}: {eq_str} = {res[i]} \u2261 0 (mod {q})")

    return x

def main():
    # Define problem parameters that satisfy the condition m > n*k
    n_param = 3
    k_param = 3
    q_param = 2**k_param # q=8
    m_param = n_param * k_param + 1 # m=10, so m>nk=9

    print("--- Generating a random problem instance ---")
    # A is a random n x m matrix with entries in Z_q
    np.random.seed(42) # For reproducibility
    A_matrix = np.random.randint(0, q_param, size=(n_param, m_param))
    print("Using matrix A:")
    print(A_matrix)

    # Solve the problem
    solve_binary_hsis(A_matrix, q_param)

if __name__ == '__main__':
    main()