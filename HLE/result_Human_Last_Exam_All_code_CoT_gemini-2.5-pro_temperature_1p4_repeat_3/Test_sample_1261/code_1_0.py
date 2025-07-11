import numpy as np

def null_space_basis_F2(A):
    """
    Finds a basis for the null space of a matrix A over F_2 using Gaussian elimination.
    """
    if A.size == 0:
        return np.array([[]])
        
    m, n = A.shape
    A = A.copy()
    
    # Forward elimination to get row echelon form
    pivot_row = 0
    pivot_cols = []
    for j in range(n): # col
        if pivot_row < m:
            i = pivot_row
            while i < m and A[i, j] == 0:
                i += 1
            
            if i < m:
                A[[pivot_row, i]] = A[[i, pivot_row]] # Swap rows
                pivot_cols.append(j)
                # Eliminate other 1s in the same column
                for i in range(m):
                    if i != pivot_row and A[i, j] == 1:
                        A[i, :] = (A[i, :] + A[pivot_row, :]) % 2
                pivot_row += 1

    # Backward substitution is not strictly needed for nullspace basis
    
    # Identify free columns
    free_cols = [j for j in range(n) if j not in pivot_cols]
    
    # Construct basis for null space
    basis = []
    for free_col in free_cols:
        x = np.zeros(n, dtype=int)
        x[free_col] = 1
        for i, pivot_col in enumerate(pivot_cols):
            x[pivot_col] = A[i, free_col]
        basis.append(x)
        
    if not basis:
        return np.zeros((n, 0), dtype=int)
        
    return np.array(basis).T

def solve_Ax_zero_mod_q(A, n, m, k):
    """
    Finds a non-zero vector x in {0,1}^m such that Ax = 0 (mod 2^k).
    """
    q = 2**k
    print(f"Searching for non-zero x in {{0,1}}^{m} such that Ax = 0 (mod {q})")
    
    # B_prev columns form a basis for solutions mod 2**(i-1)
    # Initially, any vector is a solution mod 2^0=1
    B_prev = np.identity(m, dtype=int)

    for i in range(1, k + 1):
        # M_int = A @ B_prev
        M_int = np.dot(A, B_prev)
        
        # All entries in M_int should be divisible by 2**(i-1)
        divisor = 2**(i - 1)
        if not np.all((M_int % divisor) == 0):
            print(f"Error at step {i}: Not all elements are divisible by {divisor}")
            return None

        # V is the matrix for the linear system over F_2
        V = (M_int // divisor) % 2
        
        # C's columns form a basis for the nullspace of V
        C = null_space_basis_F2(V)

        if C.shape[1] == 0:
            print(f"No non-trivial solution found at step {i}. This should not happen with the given parameters.")
            return None

        # Update the basis for the next iteration
        # B_curr = (B_prev @ C) % 2
        B_curr = np.dot(B_prev, C) % 2
        B_prev = B_curr

    # Any column of the final B is a solution. We take the first one.
    x = B_prev[:, 0]
    
    return x

if __name__ == '__main__':
    # Set parameters for the problem
    # Note: For demonstration, n and k are small.
    # m = Omega(n^k) is the condition.
    n = 2
    k = 3 # q = 8. k > 1
    # We need m > n*k = 2*3 = 6.
    # We also need m = Omega(n^k) = Omega(2^3) = Omega(8). Let's take m=10.
    m = 10
    q = 2**k

    # Check if parameters are valid for the theory
    if m <= n * k:
        print(f"Warning: m={m} is not > n*k={n*k}. The algorithm is not guaranteed to succeed.")
    
    # Generate a random matrix A from Z_q^{n x m}
    A = np.random.randint(0, q, size=(n, m))
    
    print("Matrix A (mod q):")
    print(A)
    print("-" * 20)
    
    # Find the solution vector x
    x = solve_Ax_zero_mod_q(A, n, m, k)
    
    if x is not None:
        print("-" * 20)
        print(f"Found a non-zero solution vector x in {{0,1}}^{m}:")
        print(x)
        print("-" * 20)

        # Verification
        print("Verifying the solution Ax = 0 (mod q):")
        Ax = np.dot(A, x)
        Ax_mod_q = Ax % q
        
        print(f"Ax = {Ax}")
        print(f"Ax mod {q} = {Ax_mod_q}")

        is_correct = np.all(Ax_mod_q == 0)
        print(f"Solution is correct: {is_correct}")

        if is_correct:
            print("\nPrinting the final equations:")
            for i in range(n):
                equation_parts = []
                for j in range(m):
                    equation_parts.append(f"{A[i, j]}*{x[j]}")
                equation_str = " + ".join(equation_parts)
                result_val = Ax[i]
                print(f"Row {i}: {equation_str} = {result_val} \u2261 0 (mod {q})")
