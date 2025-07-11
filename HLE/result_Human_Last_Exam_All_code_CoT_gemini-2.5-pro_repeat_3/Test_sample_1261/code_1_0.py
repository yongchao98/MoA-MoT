import numpy as np

def find_null_space_basis_gf2(matrix):
    """
    Finds a basis for the null space of a matrix over GF(2) using Gaussian elimination.
    
    Args:
        matrix (np.ndarray): The input matrix with integer entries.
        
    Returns:
        list of np.ndarray: A list of basis vectors for the null space.
    """
    mat = matrix.copy() % 2
    n_rows, n_cols = mat.shape
    
    pivot_row = 0
    pivot_cols = []
    
    # Forward elimination to bring matrix to row echelon form
    for j in range(n_cols):
        if pivot_row < n_rows:
            i = pivot_row
            while i < n_rows and mat[i, j] == 0:
                i += 1
            
            if i < n_rows:
                mat[[pivot_row, i]] = mat[[i, pivot_row]] # Swap rows to bring pivot to front
                
                # Eliminate other 1s in the same column
                for other_row in range(n_rows):
                    if other_row != pivot_row and mat[other_row, j] == 1:
                        mat[other_row, :] = (mat[other_row, :] + mat[pivot_row, :]) % 2
                
                pivot_cols.append(j)
                pivot_row += 1
        else:
            break

    # Identify free variables and construct basis vectors
    basis = []
    free_vars_indices = [j for j in range(n_cols) if j not in pivot_cols]
    
    for free_var_idx in free_vars_indices:
        solution_vector = np.zeros(n_cols, dtype=int)
        solution_vector[free_var_idx] = 1
        
        # Back-substitute to find values for pivot variables
        for i, pivot_col in enumerate(pivot_cols):
            if mat[i, free_var_idx] == 1:
                solution_vector[pivot_col] = 1
            
        basis.append(solution_vector)
        
    return basis

def solve_binary_null_space(A, q):
    """
    Finds a non-zero binary vector x such that Ax = 0 (mod q), using a lifting algorithm.
    
    Args:
        A (np.ndarray): The n x m matrix.
        q (int): The modulus, a power of 2 greater than 2.
        
    Returns:
        np.ndarray or None: The solution vector x, or None if not found.
    """
    n, m = A.shape
    k = q.bit_length() - 1

    # Step 1: Find basis for solutions mod 2
    current_basis = find_null_space_basis_gf2(A)
    
    if not current_basis:
        print("No non-trivial solution mod 2 exists.")
        return None

    modulus = 2
    
    # Step 2: Iteratively lift the solution from mod 2^j to mod 2^(j+1)
    for _ in range(1, k):
        modulus *= 2
        
        if not current_basis:
            print(f"No solution space to lift from for modulus {modulus}")
            return None

        # Build the matrix for the next system of equations
        d = len(current_basis)
        V = np.zeros((n, d), dtype=int)
        
        for i, b in enumerate(current_basis):
            Ab = A @ b
            v = (Ab // (modulus // 2)) % 2
            V[:, i] = v
            
        # Find the null space for the coefficients of the linear combination
        coeff_basis = find_null_space_basis_gf2(V)
        
        if not coeff_basis:
            # This is unexpected given the problem's constraints (m >> kn)
            print(f"Lifting failed at modulus {modulus}. No non-trivial solution found.")
            return None

        # Construct the new basis for the lifted solutions
        next_basis = []
        for c in coeff_basis:
            new_b = np.zeros(m, dtype=int)
            # Combine old basis vectors according to the new coefficient vector c
            for i_vec, coeff_val in enumerate(c):
                if coeff_val == 1:
                    new_b = (new_b + current_basis[i_vec])
            new_b %= 2
            if np.any(new_b): # Ensure the new basis vector is non-zero
                 next_basis.append(new_b)
        
        current_basis = next_basis

    return current_basis[0] if current_basis else None

# --- Main execution demonstrating the algorithm ---
# Setup parameters according to the problem description
n = 3
k = 2  # k > 1
q = 2**k # q = 4

# m must be Omega(n^k) and poly(n).
# n^k = 3^2 = 9. Let's choose m > 9.
m = 10 
print(f"Parameters: n={n}, m={m}, k={k}, q={q}")
print(f"Constraint check: m={m}, n^k={n**k}. (m=Omega(n^k) is satisfied)")

# Generate a random matrix A from Z_q^{n x m}
np.random.seed(0) # for reproducibility
A = np.random.randint(0, q, size=(n, m))

# Find the solution vector x
x = solve_binary_null_space(A, q)

if x is not None:
    print("\nFound a non-zero solution vector x.")
    result_vector = (A @ x) % q
    
    print(f"\nMatrix A (sample):\n{A[:,:8]}\n...")
    print(f"\nSolution x:\n{x}")
    print(f"\nModulus q: {q}")
    
    print("\nFinal Equation Verification: A * x (mod q)")
    print(f"(A @ x) % {q} = {result_vector}")
    
    if np.all(result_vector == 0):
        print("\nSuccess! The equation Ax = 0 (mod q) holds.")
    else:
        print("\nFailure! The equation does not hold.")
else:
    print("\nNo solution was found (unexpected).")
