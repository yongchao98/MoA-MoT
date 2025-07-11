import numpy as np

def find_f2_null_space_basis(matrix):
    """
    Finds a basis for the null space of a matrix over F_2 using Gaussian elimination.
    The basis vectors are returned as columns of the resulting matrix.
    """
    M = np.copy(matrix).astype(int) % 2
    n_rows, n_cols = M.shape
    pivot_row = 0
    pivot_cols = []

    # Forward elimination to get row echelon form
    for j in range(n_cols):
        if pivot_row >= n_rows:
            break
        i = pivot_row
        while i < n_rows and M[i, j] == 0:
            i += 1
        if i < n_rows:
            M[[pivot_row, i]] = M[[i, pivot_row]]
            pivot_cols.append(j)
            for r in range(pivot_row + 1, n_rows):
                if M[r, j] == 1:
                    M[r, :] = (M[r, :] + M[pivot_row, :]) % 2
            pivot_row += 1

    # Backward substitution to get reduced row echelon form (RREF)
    for i in range(len(pivot_cols) - 1, -1, -1):
        pivot_col = pivot_cols[i]
        for r in range(i):
            if M[r, pivot_col] == 1:
                M[r, :] = (M[r, :] + M[i, :]) % 2
    
    RREF = M
    free_cols = [j for j in range(n_cols) if j not in pivot_cols]
    basis = []
    for free_col in free_cols:
        x = np.zeros(n_cols, dtype=int)
        x[free_col] = 1
        for i, pivot_col in enumerate(pivot_cols):
            if RREF[i, free_col] == 1:
                x[pivot_col] = 1
        basis.append(x)
    
    return np.array(basis).T if basis else np.zeros((n_cols, 0), dtype=int)

def demonstrate_lifting_algorithm():
    """
    Solves Ax = 0 (mod q) for a non-zero x in {0,1}^m on a toy example
    to demonstrate the deterministic lifting algorithm.
    """
    # Setup a toy problem: n=1, k=2 (q=4), m=3.
    # The condition m > n*k holds (3 > 1*2).
    n, k_val, m, q = 1, 2, 3, 4
    # Let A be a 1x3 matrix.
    A = np.array([[1, 2, 3]])

    print("Demonstrating the deterministic lifting algorithm.")
    print(f"Goal: Find non-zero binary vector x such that Ax = 0 (mod {q})")
    print(f"Matrix A = {A}")
    print("-" * 30)

    # Step 1: Solve Ax = 0 (mod 2). This gives a basis for L_1.
    print("Step 1: Solve Ax = 0 (mod 2)")
    A_mod2 = A % 2
    # The system is x_1 + x_3 = 0 (mod 2)
    basis_L1 = find_f2_null_space_basis(A_mod2)
    print(f"Basis for solutions mod 2 (L1) are columns of:\n{basis_L1}")
    
    # Step 2: Lift solutions from mod 2 to mod 4.
    print("\nStep 2: Lift solutions to mod 4")
    # For each basis vector b in L1, A*b = 0 (mod 2), so A*b = 2*v for some integer vector v.
    # We want to find a linear combination x = L1*c such that Ax = 0 (mod 4).
    # This is equivalent to finding c such that [(A*L1)/2]*c = 0 (mod 2).
    V = (A @ basis_L1) // 2
    print(f"Constructing new system matrix V = (A @ L1) / 2:\n{V}")
    V_mod2 = V % 2
    
    # Solve Vc = 0 (mod 2) for coefficients c
    print("Solving Vc = 0 (mod 2)...")
    c_basis = find_f2_null_space_basis(V_mod2)
    
    # Any non-zero column in c_basis gives a valid coefficient vector.
    c = c_basis[:, 0]
    print(f"Found coefficients c = {c}")

    # The final solution x is the combination of L1 basis vectors.
    # This procedure is guaranteed to find a {0,1} solution given the problem constraints.
    # In this case, the resulting vector is already in {0,1}^m.
    x_final = (basis_L1 @ c)

    print("\n--- Final Equation ---")
    print(f"A = {A}")
    print(f"x = {x_final}")
    result_vec = (A @ x_final)
    result_mod_q = result_vec % q
    print(f"A x = {A} @ {x_final} = {result_vec[0]}")
    print(f"A x (mod {q}) = {result_mod_q[0]}")

    is_correct = (result_mod_q[0] == 0)
    print(f"\nVerification: Is Ax=0 (mod {q})? {is_correct}")

demonstrate_lifting_algorithm()