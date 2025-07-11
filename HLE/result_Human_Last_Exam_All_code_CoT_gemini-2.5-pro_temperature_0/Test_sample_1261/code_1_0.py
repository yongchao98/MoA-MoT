import numpy as np

def find_null_space_basis_z2(M):
    """
    Finds a basis for the null space of a matrix M over Z_2 using Gaussian elimination.
    """
    A = np.copy(M).astype(int)
    n_rows, n_cols = A.shape
    pivot_row = 0
    pivot_cols = []

    for j in range(n_cols):
        if pivot_row < n_rows:
            pivot = np.where(A[pivot_row:, j] == 1)[0]
            if pivot.size > 0:
                pivot_actual_row = pivot_row + pivot[0]
                A[[pivot_row, pivot_actual_row]] = A[[pivot_actual_row, pivot_row]]
                
                for i in range(n_rows):
                    if i != pivot_row and A[i, j] == 1:
                        A[i, :] = (A[i, :] + A[pivot_row, :]) % 2
                
                pivot_cols.append(j)
                pivot_row += 1

    rank = len(pivot_cols)
    free_cols = [j for j in range(n_cols) if j not in pivot_cols]
    
    basis = []
    for free_col_idx in free_cols:
        b = np.zeros(n_cols, dtype=int)
        b[free_col_idx] = 1
        for i, pivot_col_idx in reversed(list(enumerate(pivot_cols))):
            val = A[i, free_col_idx]
            if val == 1:
                b[pivot_col_idx] = 1
        basis.append(b)
        
    if not basis:
        return np.zeros((n_cols, 0), dtype=int)

    return np.array(basis).T

def solve_mod_power_of_2(A, q, n, m, k):
    """
    Finds a non-zero x in {0,1}^m such that Ax = 0 (mod q), where q=2^k.
    """
    # Step 1: Initialization for mod 2
    A_mod_2 = A % 2
    S = find_null_space_basis_z2(A_mod_2)

    if S.shape[1] == 0:
        # This should not happen given the problem constraints m > n
        return None

    # Step 2: Iterative Lifting
    for i in range(1, k):
        power_of_2 = 2**i
        
        # A is int, S is int {0,1}. Compute A@S over integers.
        AS = A @ S
        
        # Check if divisible, then perform integer division
        if not np.all((AS % power_of_2) == 0):
             raise Exception(f"Lifting failed at step i={i}. Not all elements of A@S are divisible by {power_of_2}")

        T = (AS // power_of_2) % 2
        
        K = find_null_space_basis_z2(T)
        
        if K.shape[1] == 0:
            # This should not happen given m = Omega(n^k)
            return None
            
        S = (S @ K) % 2

    # Step 3: Final Solution
    # Any non-zero column of the final S is a solution
    if S.shape[1] > 0:
        # Find first non-zero column
        for j in range(S.shape[1]):
            x = S[:, j]
            if np.any(x):
                return x
    return None

if __name__ == '__main__':
    # Set parameters according to the problem statement
    # k > 1, m = Omega(n^k), m in poly(n)
    n = 3
    k = 2
    # Let's choose m = n^k + 1 to satisfy the constraints
    m = n**k + 1
    q = 2**k

    print(f"Parameters: n={n}, m={m}, k={k}, q={q}\n")

    # Generate a random matrix A from Z_q^{n x m}
    # For reproducibility, we can set a seed
    np.random.seed(0)
    A = np.random.randint(0, q, size=(n, m))

    # Find the solution
    x = solve_mod_power_of_2(A, q, n, m, k)

    if x is not None:
        print("Found a non-zero solution x.")
        # Verification
        res = (A @ x) % q
        
        print("The equation is Ax = 0 (mod q), where:")
        print("\nA =")
        print(A)
        print("\nx =")
        print(x)
        print(f"\nAx (mod {q}) =")
        print(res)

        if np.all(res == 0):
            print("\nVerification successful: Ax is indeed 0 (mod q).")
        else:
            print("\nVerification failed.")
    else:
        print("No non-zero solution was found.")