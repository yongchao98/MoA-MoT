import numpy as np

def find_null_space_Z2(A):
    """
    Finds a basis for the null space of a matrix A over Z_2 using Gaussian elimination.
    """
    n, m = A.shape
    M = np.copy(A).astype(int)
    
    pivot_row = 0
    pivot_cols = []
    
    # Forward elimination to get row echelon form
    for col in range(m):
        if pivot_row < n:
            pivot = np.where(M[pivot_row:, col] == 1)[0]
            if pivot.size > 0:
                pivot_actual_row = pivot[0] + pivot_row
                # Swap rows to bring pivot to the top of the submatrix
                M[[pivot_row, pivot_actual_row]] = M[[pivot_actual_row, pivot_row]]
                # Eliminate other 1s in the same column
                for row in range(n):
                    if row != pivot_row and M[row, col] == 1:
                        M[row] = (M[row] + M[pivot_row]) % 2
                pivot_cols.append(col)
                pivot_row += 1

    # Identify free columns
    free_cols = [j for j in range(m) if j not in pivot_cols]
    
    # Back substitution to find basis vectors for the null space
    null_space_basis = []
    for free_col in free_cols:
        basis_vector = np.zeros(m, dtype=int)
        basis_vector[free_col] = 1
        for i, pivot_col in reversed(list(enumerate(pivot_cols))):
            val = M[i, free_col]
            if val == 1:
                basis_vector[pivot_col] = 1
        
        # In Z2, we need to solve M_pivot * x_pivot = M_free * x_free
        # The logic is slightly different for reduced row echelon form
    
    # Simplified approach: create reduced row echelon form first
    M = np.copy(A).astype(int)
    pivot_row = 0
    pivot_cols = []
    col_to_pivot_row = {}

    for col in range(m):
        if pivot_row < n:
            i = pivot_row
            while i < n and M[i, col] == 0:
                i += 1
            if i < n: # Found pivot
                M[[pivot_row, i]] = M[[i, pivot_row]]
                for j in range(n):
                    if j != pivot_row and M[j, col] == 1:
                        M[j] = (M[j] + M[pivot_row]) % 2
                pivot_cols.append(col)
                col_to_pivot_row[col] = pivot_row
                pivot_row += 1

    free_cols = [j for j in range(m) if j not in pivot_cols]
    
    basis = []
    for free_col in free_cols:
        v = np.zeros(m, dtype=int)
        v[free_col] = 1
        for pivot_col, p_row in col_to_pivot_row.items():
            if M[p_row, free_col] == 1:
                v[pivot_col] = 1
        basis.append(v)
        
    return basis

def solve_sis_instance(n, k, m):
    """
    Demonstrates the lifting algorithm for a random instance.
    """
    if m <= n * k:
        print(f"Error: m must be greater than n*k. Got m={m}, n*k={n*k}.")
        return

    q = 2**k
    print(f"Solving for n={n}, k={k} (q={q}), m={m}")
    print("-" * 30)

    # 1. Generate a random matrix A
    A = np.random.randint(0, q, size=(n, m))
    print("Random matrix A (first 5 columns):\n", A[:,:5])
    print("-" * 30)

    # 2. Base Case: Find null space basis for A mod 2
    A_mod_2 = A % 2
    current_basis = find_null_space_Z2(A_mod_2)
    
    if not current_basis:
        print("Could not find a non-trivial null space mod 2. This is unexpected.")
        return

    # 3. Lifting Loop
    for j in range(1, k):
        p_j = 2**j
        s = len(current_basis)
        
        # Build the matrix M for the next linear system
        M = np.zeros((n, s), dtype=int)
        for i, v in enumerate(current_basis):
            Av = A @ v
            M[:, i] = (Av // p_j) % 2
        
        # Find the null space of M to find the coefficients for the new basis
        kernel_M = find_null_space_Z2(M)
        
        if not kernel_M:
            print(f"Lifting failed at step j={j}. Could not find a non-trivial kernel.")
            return

        # Construct the new basis for solutions mod 2^(j+1)
        new_basis = []
        for c in kernel_M:
            v_new = np.zeros(m, dtype=int)
            for i_c, c_val in enumerate(c):
                if c_val == 1:
                    v_new = (v_new + current_basis[i_c]) % 2
            new_basis.append(v_new)
        current_basis = new_basis

    # 4. Final Solution
    # Any vector from the final basis is a solution
    x = current_basis[0]

    # 5. Verification
    result_vector = (A @ x) % q
    is_solution = np.all(result_vector == 0)

    print(f"Found a non-zero binary vector x:\n{x}")
    print("\nVerification:")
    # We print the equation Ax = 0 mod q
    for i in range(n):
      line = " + ".join([f"{A[i, j]}*{x[j]}" for j in range(m)])
      print(f"({line}) mod {q} = {result_vector[i]}")
    
    print(f"\nIs Ax = 0 (mod q)? {is_solution}")


if __name__ == '__main__':
    # Problem parameters
    n = 2  # Number of equations
    k = 3  # Power for q = 2^k
    # m must be > n*k to guarantee a solution at each lifting step
    m = n * k + 1 
    
    solve_sis_instance(n, k, m)