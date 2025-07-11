import numpy as np

def rref_gf2(matrix):
    """
    Computes the Row-Reduced Echelon Form of a matrix over GF(2).
    Returns the RREF matrix and a list of pivot column indices.
    """
    mat = matrix.copy().astype(int)
    n_rows, n_cols = mat.shape
    pivot_row = 0
    pivot_cols = []
    for j in range(n_cols): # Iterate through columns
        if pivot_row < n_rows:
            pivot_found = False
            for i in range(pivot_row, n_rows):
                if mat[i, j] == 1:
                    mat[[pivot_row, i]] = mat[[i, pivot_row]] # Swap rows
                    pivot_found = True
                    break
            
            if pivot_found:
                pivot_cols.append(j)
                for i in range(n_rows):
                    if i != pivot_row and mat[i, j] == 1:
                        mat[i, :] = (mat[i, :] + mat[pivot_row, :]) % 2
                pivot_row += 1
    return mat, pivot_cols

def null_space_gf2(matrix):
    """
    Computes a basis for the null space of a matrix over GF(2).
    Returns a matrix whose columns form the basis.
    """
    n, m = matrix.shape
    rref, pivot_cols = rref_gf2(matrix)
    free_cols = [j for j in range(m) if j not in pivot_cols]
    
    d = len(free_cols)
    if d == 0:
        return np.zeros((m, 0), dtype=int)
        
    basis = np.zeros((m, d), dtype=int)
    
    for i, free_col_idx in enumerate(free_cols):
        basis[free_col_idx, i] = 1
        for row_idx, pivot_col_idx in enumerate(pivot_cols):
            basis[pivot_col_idx, i] = rref[row_idx, free_col_idx]
            
    return basis

def find_solution_recursive(A, k):
    """
    Recursively finds a non-zero binary solution to Ax = 0 (mod 2^k).
    """
    n, m = A.shape
    if k == 0:
        # Base case: mod 1. Any non-zero vector is a solution.
        # Return the simplest one.
        x = np.zeros(m, dtype=int)
        if m > 0:
            x[0] = 1
        return x

    # Find basis for the null space of A mod 2
    V = null_space_gf2(A % 2)
    
    # Dimension of the null space
    d = V.shape[1]
    if d == 0:
        # This case should not be reached given the problem's constraints
        return None 

    # Lifting step
    # We are looking for a solution x = V*c
    # A * (V*c) = 0 (mod 2^k)  => (A*V)*c = 0 (mod 2^k)
    # Let A_next = A*V. We know A_next is divisible by 2.
    # (A_next/2)*c = 0 (mod 2^{k-1})
    A_next = (A @ V) // 2
    
    # Recursively solve for c
    c = find_solution_recursive(A_next, k - 1)
    
    if c is None:
        return None
        
    # Reconstruct the solution x from c.
    # The multiplication is over GF(2).
    x = (V @ c) % 2
    
    return x

def solve_sis_power_of_2(n, k_power):
    """
    Generates a problem instance and solves it.
    """
    if not (isinstance(n, int) and n > 0):
        print("Error: n must be a positive integer.")
        return
    if not (isinstance(k_power, int) and k_power > 1):
        print("Error: k must be an integer > 1.")
        return

    # m must be poly(n) and Omega(n^k)
    # For demonstration, we set m > n*k to ensure a solution exists.
    # A simple choice satisfying the constraints is m = n^k + 1 for larger n,
    # or just m = n * k_power + 1 for a working example.
    m = n * k_power + 5 # Add a small margin to make it very likely non-trivial
    q = 2**k_power

    print(f"Problem setup:")
    print(f"n = {n}")
    print(f"k = {k_power}, so q = 2^k = {q}")
    print(f"m = {m} (we need m > n*k = {n*k_power})")
    print("-" * 20)

    # Generate a random matrix A
    A = np.random.randint(0, q, size=(n, m))
    
    print("Finding a non-zero binary vector x such that Ax = 0 (mod q)...")
    
    # Find the solution
    x = find_solution_recursive(A, k_power)

    if x is not None:
        print("\nSolution found!")
        # Ax = b (mod q)
        b = (A @ x) % q
        
        print(f"\nMatrix A (dims {n}x{m}):\n{A}")
        print(f"\nSolution x (binary, non-zero, dims {m}):\n{x}")
        print(f"\nVerification: Ax (mod {q})")
        
        # Print the equation A*x = b (mod q)
        # For readability, only show parts of the equation if matrices are large
        rows_to_show = min(n, 10)
        
        for i in range(rows_to_show):
            row_str = " + ".join([f"{A[i, j]}*{x[j]}" for j in range(min(m,10))])
            if m > 10:
                row_str += " + ..."
            print(f"  (Row {i})  {row_str} = {b[i]} (mod {q})")
        
        if np.all(b == 0):
            print("\nSuccess: Ax is the zero vector modulo q.")
        else:
            print("\nFailure: Ax is NOT the zero vector modulo q.")
    else:
        print("No solution was found.")

if __name__ == '__main__':
    # Example Parameters (feel free to change)
    # n and k should be small for the code to run fast and be readable.
    n_val = 2
    k_val = 3
    solve_sis_power_of_2(n_val, k_val)
