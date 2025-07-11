import numpy as np

def find_nullspace_basis_mod2(mat):
    """
    Finds a basis for the null space of a matrix over Z_2 using Gaussian elimination.
    """
    M = np.copy(mat) % 2
    n, m = M.shape
    pivot_row = 0
    pivot_cols = []
    
    # Forward elimination to get row echelon form
    for col in range(m):
        if pivot_row < n:
            pivot = np.where(M[pivot_row:, col] == 1)[0]
            if pivot.size > 0:
                pivot_idx = pivot[0] + pivot_row
                M[[pivot_row, pivot_idx]] = M[[pivot_idx, pivot_row]] # Swap rows
                pivot_cols.append(col)
                for i in range(n):
                    if i != pivot_row and M[i, col] == 1:
                        M[i, :] = (M[i, :] + M[pivot_row, :]) % 2
                pivot_row += 1
    
    # Identify free columns
    free_cols = [i for i in range(m) if i not in pivot_cols]
    
    # Back substitution to find basis vectors
    basis = []
    for free_col in free_cols:
        b = np.zeros(m, dtype=int)
        b[free_col] = 1
        for i, pivot_col in reversed(list(enumerate(pivot_cols))):
            val = np.dot(M[i, free_cols], b[free_cols]) % 2
            if val == 1:
                 b[pivot_col] = 1
        basis.append(b)
        
    if not basis: # Should not happen if m > n
        return np.array([]).reshape(m, 0)
        
    return np.array(basis).T

def solve_binary_kernel(A, q):
    """
    Solves Ax = 0 (mod q) for a non-zero binary vector x.
    q must be a power of 2.
    """
    if q == 1:
        # Trivial case, any non-zero x works
        x = np.zeros(A.shape[1], dtype=int)
        x[0] = 1
        return x
    
    if int(np.log2(q)) != np.log2(q):
        raise ValueError("q must be a power of 2.")

    # Step 1: Find basis for the nullspace of A mod 2
    basis_mod2 = find_nullspace_basis_mod2(A)
    
    if basis_mod2.shape[1] == 0:
        # No non-trivial solution mod 2, so no solution for mod q either
        return None

    # If q=2, we are done. Return the first basis vector.
    if q == 2:
        return basis_mod2[:, 0]
        
    # Step 2: Form the new matrix for the smaller modulus problem
    # Note: Use integer arithmetic, not modulo q
    V = A @ basis_mod2 
    
    # Check if all elements of V are divisible by 2
    if not np.all((V % 2) == 0):
       # This should not happen if basis_mod2 is correct
       raise RuntimeError("Error in lifting step: matrix V not divisible by 2.")

    Y = V // 2
    
    # Step 3: Recursively solve for the smaller modulus
    c = solve_binary_kernel(Y, q // 2)

    if c is None:
        return None

    # Step 4: Lift the solution back
    x = (basis_mod2 @ c) % 2
    return x

def main():
    """
    Main function to run a demonstration of the algorithm.
    """
    # Problem parameters satisfying the constraints
    n = 2
    k = 3
    # Let's set m = n^k to satisfy m = Omega(n^k) and m in poly(n)
    m = n**k 
    q = 2**k

    print(f"Solving Ax = 0 (mod q) for a non-zero binary vector x.")
    print(f"Parameters: n={n}, m={m}, k={k}, q={q}\n")

    # Generate a random matrix A from Z_q
    A = np.random.randint(0, q, size=(n, m))
    
    print("Generated random matrix A (first 10 columns):")
    print(A[:,:10])
    print("-" * 30)

    # Solve the problem
    x = solve_binary_kernel(A, q)

    if x is not None:
        print("Found a non-zero binary solution vector x:")
        print(x)
        print("-" * 30)
        
        # Verification
        print("Verifying the solution: Ax (mod q)")
        result = (A @ x) % q
        print(result)
        
        if np.all(result == 0):
            print("\nVerification successful: Ax is the zero vector modulo q.")
        else:
            print("\nVerification failed.")
            
        print("\nFinal Equation Form (conceptual):")
        print("A * x = 0 (mod q)")
        print(f"Where A is a {n}x{m} matrix, x is a {m}-dim binary vector, and q={q}.")
        
    else:
        print("Could not find a solution.")

if __name__ == "__main__":
    main()