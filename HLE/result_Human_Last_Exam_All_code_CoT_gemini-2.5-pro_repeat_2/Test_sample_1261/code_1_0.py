import numpy as np

def find_null_space_F2(A):
    """
    Computes a basis for the null space (kernel) of a matrix A over the field F_2
    using Gaussian elimination (Reduced Row Echelon Form).
    Returns a matrix whose columns form the basis.
    """
    # Create a copy to avoid modifying the original matrix
    A_rref = np.copy(A).astype(int)
    m, n = A_rref.shape
    pivot_row = 0
    pivot_cols = []
    
    # Compute Reduced Row Echelon Form (RREF)
    for j in range(n):  # Iterate through columns
        if pivot_row >= m:
            break
        
        i = pivot_row
        while i < m and A_rref[i, j] == 0:
            i += 1
        
        if i < m:  # Found a pivot in row i
            pivot_cols.append(j)
            A_rref[[pivot_row, i]] = A_rref[[i, pivot_row]]  # Swap rows
            
            # Eliminate other 1s in the current pivot column
            for i_ in range(m):
                if i_ != pivot_row and A_rref[i_, j] == 1:
                    A_rref[i_, :] = (A_rref[i_, :] + A_rref[pivot_row, :]) % 2
            pivot_row += 1

    # Extract basis vectors from RREF
    basis = []
    free_cols = [j for j in range(n) if j not in pivot_cols]
    
    for free_col in free_cols:
        vec = np.zeros(n, dtype=int)
        vec[free_col] = 1
        for i, pivot_col in enumerate(pivot_cols):
            vec[pivot_col] = A_rref[i, free_col]
        basis.append(vec)

    if not basis:
        return np.array([]).reshape(n, 0)
        
    return np.array(basis).T

def solve_binary_kernel(A, q, k):
    """
    Finds a non-zero x in {0,1}^m such that Ax = 0 (mod q), where q=2^k.
    This function implements the deterministic lifting algorithm.
    """
    n, m = A.shape
    
    # B_j is the basis matrix for solutions mod 2^j. It's an m x d_j matrix.
    # Start with B_0 = I_m, the basis for all m-dimensional binary vectors.
    B_j = np.eye(m, dtype=int)

    for j in range(1, k + 1):
        # We are lifting solutions from mod 2^(j-1) to mod 2^j.
        
        # We need to solve (A * B_{j-1}) * c = 0 (mod 2^j).
        # We know from the previous step that (A * B_{j-1}) is divisible by 2^(j-1).
        M_j_numerator = np.dot(A, B_j)
        
        # Integer division is guaranteed to be exact by the algorithm's logic.
        M_j = M_j_numerator // (2**(j - 1))
        
        # The congruence becomes M_j * c = 0 (mod 2).
        M_j_mod2 = M_j % 2
        
        # Find the null space (kernel) of M_j_mod2. C_j's columns are the basis.
        C_j = find_null_space_F2(M_j_mod2)
        
        if C_j.shape[1] == 0: # Null space is trivial (only the zero vector)
            print("No non-zero solution found. This shouldn't happen with the given parameters.")
            return None
            
        # Update the overall basis for the next iteration: B_j = B_{j-1} * C_j.
        # The matrix multiplication is effectively done over F_2.
        B_j = np.dot(B_j, C_j) % 2

    # After k steps, B_k is the basis for solutions to Ax=0 (mod q).
    # Any non-zero column of B_k is a valid non-zero binary solution vector x.
    if B_j.shape[1] > 0:
        x = B_j[:, 0]
        return x
    else:
        # This case is not expected given the problem's condition m=Ω(n^k)
        return None

# --- Main execution ---
# 1. Set problem parameters
k = 3 
n = 2
# Per problem, m = Ω(n^k). For the algorithm to be guaranteed to work, we need m > n*k.
m = n * k + 2  # Example: m = 2*3 + 2 = 8
q = 2**k

# 2. Generate a random matrix A
np.random.seed(42) # for reproducibility
A = np.random.randint(0, q, size=(n, m))

print(f"Let k = {k}, so q = 2^k = {q}")
print(f"Let n = {n}, m = {m}")
print(f"\nMatrix A (randomly sampled from Z_{q}^{{{n}x{m}}}):")
print(A)

# 3. Run the algorithm to find the solution vector x
x = solve_binary_kernel(A, q, k)

# 4. Print and verify the result
if x is not None:
    print("\nFound a non-zero binary vector x:")
    print(x)
    
    # Verification
    Ax = np.dot(A, x)
    Ax_mod_q = Ax % q
    
    print("\nVerification of the equation Ax = 0 (mod q):")
    for i in range(n):
        equation_str = " + ".join([f"{A[i, j]}*{x[j]}" for j in range(m)])
        print(f"Row {i+1}: {equation_str} = {Ax[i]} \u2261 {Ax_mod_q[i]} (mod {q})")
    
    if np.all(Ax_mod_q == 0):
        print("\nSuccess: Ax = 0 (mod q) is satisfied.")
    else:
        print("\nFailure: Ax is not 0 (mod q).")
else:
    print("\nCould not find a solution.")
