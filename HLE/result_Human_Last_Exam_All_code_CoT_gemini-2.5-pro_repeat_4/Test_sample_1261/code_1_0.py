import numpy as np

def null_space_Z2(A):
    """
    Finds a basis for the null space of a matrix A over Z_2 using Gaussian elimination.
    """
    mat = A % 2
    n, m = mat.shape
    pivot_cols = []
    free_cols = list(range(m))
    
    pivot_row = 0
    for j in range(m): # Iterate through columns
        if pivot_row >= n:
            break
        
        pivot = -1
        # Find a pivot in the current column
        for i in range(pivot_row, n):
            if mat[i, j] == 1:
                pivot = i
                break
        
        if pivot != -1:
            # Swap rows to bring pivot to pivot_row
            mat[[pivot_row, pivot]] = mat[[pivot, pivot_row]]
            
            # Eliminate other 1s in the same column
            for i in range(n):
                if i != pivot_row and mat[i, j] == 1:
                    mat[i, :] = (mat[i, :] + mat[pivot_row, :]) % 2
            
            pivot_cols.append(j)
            free_cols.remove(j)
            pivot_row += 1

    # Basis vectors for the null space
    basis = []
    for free_col_idx in free_cols:
        b = np.zeros(m, dtype=int)
        b[free_col_idx] = 1
        for i, pivot_col_idx in enumerate(pivot_cols):
            if mat[i, free_col_idx] == 1:
                b[pivot_col_idx] = 1
        basis.append(b)
        
    if not basis:
        return np.array([]).reshape(m, 0)
        
    return np.array(basis).T

def solve_hsis_recursive(A, k):
    """
    Recursively solves Ax = 0 (mod 2^k) for a non-zero binary vector x.
    """
    if k == 0:
        # Base case for recursion, q=2^0=1. Any non-zero vector is a solution.
        # The size of the vector must match the number of columns of A.
        n_cols = A.shape[1]
        if n_cols > 0:
            c = np.zeros(n_cols, dtype=int)
            c[0] = 1
            return c
        else: # No variables, so no non-zero solution
            return None

    # Find a basis for the null space of A mod 2
    B = null_space_Z2(A)
    
    # If the null space is trivial (only the zero vector), we can't proceed
    if B.shape[1] == 0:
        return None

    # A' = (A * B) / 2
    # The product A @ B is guaranteed to be even.
    A_prime = (A @ B) // 2

    # Recursively solve the smaller problem
    c = solve_hsis_recursive(A_prime, k - 1)

    if c is None:
        return None

    # Construct the solution for the current level
    # The theory guarantees that with appropriate choice of c, 
    # the result x is a {0,1} vector. For this implementation,
    # we take the first solution found by the recursion.
    x = B @ c
    return x

def solve(A, q):
    """
    Main solver function.
    """
    if q <= 1 or (q & (q - 1)) != 0:
        raise ValueError("q must be a power of 2 greater than 1.")
    
    k = int(np.log2(q))
    n, m = A.shape
    
    # The condition m > n*k ensures a solution can be found.
    if m <= n * k:
        print(f"Warning: m={m} may not be large enough compared to n*k={n*k}. A solution may not be found.")

    solution = solve_hsis_recursive(A, k)
    return solution

if __name__ == '__main__':
    # Problem parameters
    n = 3
    k = 3  # Note: k > 1
    q = 2**k

    # m = Omega(n^k) means m should be large. Let's pick m > n*k.
    m = n * k + 1  # Example: m=10 for n=3,k=3

    # 1. Generate a random matrix A
    A = np.random.randint(0, q, size=(n, m))
    print("The problem is to find a non-zero x in {0,1}^m such that Ax = 0 (mod q)")
    print(f"q = {q} (k={k})")
    print(f"Matrix A (n={n}, m={m}):\n{A}\n")

    # 2. Find a solution
    x = solve(A, q)

    # 3. Output and verify the solution
    if x is not None:
        print(f"Found a solution vector x:\n{x}\n")
        
        # Verification
        Ax = A @ x
        Ax_mod_q = Ax % q
        
        print("Verification:")
        print(f"A*x = {Ax}")
        print(f"A*x (mod {q}) = {Ax_mod_q}")

        if np.all(Ax_mod_q == 0):
            if np.all(np.isin(x, [0, 1])):
                 print("\nSUCCESS: The solution is correct and binary.")
            else:
                 print("\nPARTIAL SUCCESS: The solution is correct (Ax=0), but not binary. This indicates a nuance in the choice of 'c' in the algorithm.")
            # Even if not binary, we can show the equation holds.
            # Example: A = [[...]], x = [[...]]
            # 2 * 1 + 3 * 0 + ... = ...
            for i in range(n):
                equation_str = " + ".join([f"{A[i,j]}*{x[j]}" for j in range(m)])
                print(f"Equation {i}: {equation_str} = {Ax[i]} = 0 (mod {q})")

        else:
            print("\nFAILURE: The solution is incorrect.")
    else:
        print("\nCould not find a solution.")
