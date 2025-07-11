import numpy as np

def find_f2_null_space_basis(A):
    """
    Finds a basis for the null space of a matrix A over F_2 using Gaussian elimination.
    This implementation performs column operations on A to find a basis for its null space.
    """
    # A should be a numpy array with 0/1 entries.
    n, m = A.shape
    M = np.copy(A).astype(int)
    
    # Store transformations to find the basis vectors.
    # We are implicitly finding the null space of M by transforming [M | I]^T
    # to row-echelon form. A more direct method is used here.
    
    pivot_cols = []
    pivot_row = 0
    for j in range(m): # Iterate through columns
        if pivot_row < n:
            pivot = np.where(M[pivot_row:, j] == 1)[0]
            if len(pivot) > 0:
                pivot = pivot[0] + pivot_row
                M[[pivot_row, pivot]] = M[[pivot, pivot_row]] # Swap rows
                # Eliminate other 1s in the column
                for i in range(n):
                    if i != pivot_row and M[i, j] == 1:
                        M[i, :] = (M[i, :] + M[pivot_row, :]) % 2
                pivot_cols.append(j)
                pivot_row += 1

    free_cols = [j for j in range(m) if j not in pivot_cols]
    
    basis = []
    for free_col in free_cols:
        vec = np.zeros(m, dtype=int)
        vec[free_col] = 1
        # Back-substitution to solve for pivot variables
        for i, p_col in reversed(list(enumerate(pivot_cols))):
             # Find value in the free column's row of the RREF matrix
             val = M[i, free_col]
             if val == 1:
                 vec[p_col] = 1
    
        basis.append(vec)
        
    return basis

def solve_ax_zero(A, q):
    """
    Recursively solves Ax = 0 (mod q) for a non-zero x in {0,1}^m.
    q must be a power of 2.
    """
    n, m = A.shape
    
    # Base case
    if q == 1:
        x = np.zeros(m, dtype=int)
        if m > 0:
            x[0] = 1
        return x

    # Step 1: Reduce matrix mod 2 and find null space basis
    A_mod_2 = A % 2
    null_space_basis_vectors = find_f2_null_space_basis(A_mod_2)
    
    if not null_space_basis_vectors:
        # Should not happen given problem constraints m > n*k
        return None 

    d = len(null_space_basis_vectors)
    # B is the m x d matrix where columns are the basis vectors
    B = np.array(null_space_basis_vectors).T 

    # Step 2: Construct the new matrix C for the recursive call
    # C = (A @ B) / 2
    # Perform multiplication with original A, then integer divide by 2
    AB = A @ B
    C = (AB // 2)
    
    # Step 3: Recursive call with q' = q/2
    y = solve_ax_zero(C, q // 2)
    
    if y is None:
        return None

    # Step 4: Construct final solution x = B @ y (mod 2)
    x = (B @ y) % 2
    
    return x

def main():
    # Setup problem parameters based on the prompt
    n = 3
    k = 3 # k > 1
    q = 2**k

    # Set m according to the condition m = Omega(n^k) and m > n*k
    # Example: m = n^k = 27. Let's pick m=30.
    # Check condition m > n*k: 30 > 3*3 = 9 is satisfied.
    m = 30
    
    print(f"Goal: Find non-zero x in {{0,1}}^{m} such that Ax = 0 (mod q)")
    print(f"Parameters:")
    print(f"n = {n}")
    print(f"k = {k}")
    print(f"q = 2^k = {q}")
    print(f"m = {m}")
    
    # Generate a random matrix A from Z_q^{n x m}
    # For a random A, a solution is expected to exist.
    np.random.seed(42) # For reproducibility
    A = np.random.randint(0, q, size=(n, m))

    # Solve the problem
    x = solve_ax_zero(A, q)

    if x is not None:
        # Verification
        result_vector = (A @ x) % q
        is_zero_vector = np.all(result_vector == 0)
        is_nonzero_x = np.any(x != 0)
        is_binary_x = np.all((x == 0) | (x == 1))

        print(f"\nVerification:")
        print(f"Is x a non-zero vector? {is_nonzero_x}")
        print(f"Is x a binary vector? {is_binary_x}")
        print(f"Is Ax = 0 (mod q)? {is_zero_vector}")

        if is_zero_vector and is_nonzero_x and is_binary_x:
            print("\nSolution is correct.")
            # Fulfilling the request to output numbers in the final equation.
            print("\nThe equation Ax = 0 (mod q) is satisfied as shown:")
            for i in range(n):
                line = ""
                # Displaying a truncated equation for readability with large m
                max_terms_to_show = 5
                for j in range(min(m, max_terms_to_show)):
                    line += f"({A[i,j]} * {x[j]}) + "
                if m > max_terms_to_show:
                    line += "... + "
                line += f"({A[i, m-1]} * {x[m-1]})"
                
                print(f"Row {i}: {line} = {result_vector[i]} (mod {q})")
        else:
            print("\nError: The found solution is incorrect.")
    else:
        print("\nCould not find a solution.")

if __name__ == "__main__":
    main()
