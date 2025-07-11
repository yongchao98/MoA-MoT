import numpy as np

def find_null_space_mod2(A):
    """
    Finds a basis for the null space of a binary matrix A over Z_2
    using Gaussian elimination.
    """
    rows, cols = A.shape
    A_copy = np.copy(A).astype(int)

    pivot_row = 0
    pivot_cols = []
    
    # Forward elimination
    for j in range(cols):
        if pivot_row >= rows:
            break
        i = pivot_row
        while i < rows and A_copy[i, j] == 0:
            i += 1
        
        if i < rows:
            A_copy[[pivot_row, i]] = A_copy[[i, pivot_row]]
            for k in range(rows):
                if k != pivot_row and A_copy[k, j] == 1:
                    A_copy[k, :] = (A_copy[k, :] + A_copy[pivot_row, :]) % 2
            pivot_cols.append(j)
            pivot_row += 1

    # Find free columns
    free_cols = [j for j in range(cols) if j not in pivot_cols]
    
    if not free_cols:
        return np.zeros((cols, 0), dtype=int)

    # Build null space basis
    basis = []
    for free_col_idx in free_cols:
        vec = np.zeros(cols, dtype=int)
        vec[free_col_idx] = 1
        for i, pivot_col_idx in enumerate(pivot_cols):
            if A_copy[i, free_col_idx] == 1:
                vec[pivot_col_idx] = 1
        basis.append(vec)
    
    return np.array(basis).T

def solve_binary_sis(A, q, k):
    """
    Solves Ax = 0 (mod q) for a non-zero integer vector x using a lifting algorithm.
    Under the condition m > nk, this finds a small integer solution.
    More advanced algorithms can guarantee a {0,1} solution.
    
    Args:
        A (np.array): The n x m matrix.
        q (int): The modulus, 2^k.
        k (int): The exponent.

    Returns:
        np.array: A solution vector x, or None if not found.
    """
    n, m = A.shape
    if m <= n * k:
        print(f"Algorithm requires m > n*k, but m={m} and n*k={n*k}. May fail.")

    # B's columns will span the solution space at each level of lifting.
    # We start with the identity matrix, representing all possible vectors.
    B = np.identity(m, dtype=np.int64)

    for i in range(1, k + 1):
        # At step i, B's columns are solutions mod 2**(i-1)
        # We need to find combinations of B's columns that are solutions mod 2**i
        
        power_of_2 = 1 << (i - 1)  # This is 2**(i-1)
        
        # By construction, A @ B is a matrix of multiples of power_of_2
        C = (A @ B) // power_of_2
        
        # We want to find a linear combination of columns of B, let's call it B@z,
        # such that A @ (B@z) = 0 (mod 2**i).
        # This simplifies to C@z = 0 (mod 2).
        C_mod2 = C % 2
        
        # Z's columns form a basis for the null space of C_mod2
        Z = find_null_space_mod2(C_mod2)
        
        if Z.shape[1] == 0: # If null space is trivial
            print(f"Algorithm failed at step {i}: could not find non-trivial solution.")
            return None

        # Update B. The new B's columns are solutions mod 2**i.
        B = B @ Z
    
    # Any column of the final B is a solution to Ax = 0 (mod q).
    # Take the first column as our solution vector x.
    x = B[:, [0]]
    
    return x

# --- Main execution ---
# Setup parameters satisfying m > nk and k > 1
n = 2
k = 3
q = 1 << k  # q = 2^k = 8
m = n * k + 1  # m=7, which satisfies m > nk (7 > 6)

# Generate a random matrix A
np.random.seed(42) # for reproducibility
A = np.random.randint(0, q, size=(n, m))

print("This script solves Ax = 0 (mod q) for a non-zero vector x.")
print("It uses a deterministic lifting algorithm that runs in polynomial time in n, m, k.")
print("-" * 20)
print(f"Parameters: n={n}, m={m}, k={k}, q={q}")
print(f"Condition m > nk is satisfied: {m} > {n*k}")
print("-" * 20)

# Find the solution
x = solve_binary_sis(A, q, k)

if x is not None:
    # Verify the solution
    result = (A @ x) % q
    is_solution = np.all(result == 0)

    print("Input Matrix A:\n", A)
    print("\nFound solution vector x:\n", x)
    
    # Output the final equation
    print("\nFinal equation check: Ax (mod q) = 0")
    
    a_str_list = []
    for r in range(A.shape[0]):
      row_str = " ".join([f"{val: >2}" for val in A[r, :]])
      a_str_list.append(f"[{row_str} ]")
    a_matrix_str = "\n".join(a_str_list)

    x_str_list = [f"[{val[0]: >2}]" for val in x]
    x_matrix_str = "\n".join(x_str_list)

    res_str_list = [f"[{val[0]: >2}]" for val in result]
    res_matrix_str = "\n".join(res_str_list)

    zero_str_list = [f"[{0: >2}]" for _ in range(n)]
    zero_matrix_str = "\n".join(zero_str_list)

    print(f"\n{a_matrix_str}   *   {x_matrix_str}   (mod {q})   =   {res_matrix_str}   (which should be {zero_matrix_str})")

    print(f"\nVerification successful: {is_solution}")

else:
    print("Could not find a solution.")
