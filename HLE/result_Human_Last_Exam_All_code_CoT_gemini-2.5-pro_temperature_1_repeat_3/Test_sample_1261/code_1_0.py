import numpy as np

def find_null_space_F2(A):
    """
    Finds a basis for the null space of matrix A over F_2 using Gaussian elimination.
    """
    mat = np.copy(A).astype(int) % 2
    n, m = mat.shape
    pivot_row = 0
    pivot_cols = []
    
    # Forward elimination
    for j in range(m):
        if pivot_row < n:
            pivot = np.where(mat[pivot_row:, j] == 1)[0]
            if pivot.size > 0:
                pivot = pivot[0] + pivot_row
                # Swap rows
                mat[[pivot_row, pivot], :] = mat[[pivot, pivot_row], :]
                pivot_cols.append(j)
                # Eliminate other 1s in the column
                for i in range(n):
                    if i != pivot_row and mat[i, j] == 1:
                        mat[i, :] = (mat[i, :] + mat[pivot_row, :]) % 2
                pivot_row += 1
    
    # Find free columns
    free_cols = [j for j in range(m) if j not in pivot_cols]
    
    # Basis vectors for null space
    basis = []
    for free_col in free_cols:
        vec = np.zeros(m, dtype=int)
        vec[free_col] = 1
        for i, pivot_col in enumerate(pivot_cols):
            if mat[i, free_col] == 1:
                vec[pivot_col] = 1
        basis.append(vec)
        
    if not basis: # This case should not happen given the problem constraints
        return np.zeros((m, 0), dtype=int)

    return np.array(basis).T

def solve_binary_homogeneous(A, q):
    """
    Recursively finds a non-zero x in {0,1}^m such that Ax = 0 (mod q).
    q must be a power of 2.
    """
    if q == 1:
        # Any non-zero binary vector is a solution.
        sol = np.zeros(A.shape[1], dtype=int)
        sol[0] = 1
        return sol

    if not (q > 0 and (q & (q - 1) == 0)):
        raise ValueError("q must be a power of 2.")

    # Base case k=1 (q=2)
    if q == 2:
        null_space_basis = find_null_space_F2(A)
        if null_space_basis.shape[1] == 0:
             # Should not happen given m > n
            return np.zeros(A.shape[1], dtype=int)
        # Return the first basis vector
        return null_space_basis[:, 0]

    # Recursive step
    # Find basis B for Null(A mod 2)
    B = find_null_space_F2(A)
    if B.shape[1] == 0:
        return np.zeros(A.shape[1], dtype=int)
    
    # Calculate A_next = (A @ B) / 2
    # The product A @ B has all entries even
    A_next = (A @ B) // 2
    
    # Recursively solve for z
    z = solve_binary_homogeneous(A_next, q // 2)
    
    if np.all(z == 0):
        # Sub-problem failed to find a non-zero solution
        return z
        
    # Our solution is x = B @ z
    x = B @ z
    return x

def main():
    """
    Main function to run the algorithm on an example.
    """
    # Let q= 2^k for some integer k >1. Let m= Ω(n^k), m ∈ poly(n).
    # Example: n=1, k=2 => q=4. m = Ω(1^2=1), so m can be > 1. Let's take m=3.
    n = 1
    m = 3
    k = 2
    q = 2**k

    # A is a matrix sampled uniformly random from Z_q^{n×m}
    # For demonstration, we use a fixed matrix.
    # A = [a1, a2, a3] where a_i are from Z_4.
    # E.g., A = [[1, 2, 3]]
    A = np.array([[1, 2, 3]])
    print(f"Input Matrix A:\n{A}")
    print(f"Modulus q: {q}")
    print("-" * 20)

    # Find a non-zero vector x in {0,1}^m such that Ax = 0 (mod q)
    x = solve_binary_homogeneous(A, q)

    print(f"Found solution vector x:\n{x}")

    # Verification
    res = (A @ x) % q
    is_binary = np.all((x == 0) | (x == 1))
    is_non_zero = np.any(x != 0)

    print(f"\nVerification:")
    print(f"A @ x (mod q) = {res}")
    print(f"Is x a binary vector? {is_binary}")
    print(f"Is x a non-zero vector? {is_non_zero}")

    if is_binary and is_non_zero and np.all(res == 0):
        print("\nSolution is correct.")
        print("Final equation:")
        equation_parts = []
        for i in range(m):
            if x[i] != 0:
                equation_parts.append(f"{x[i]} * {A[0, i]}")
        print(" + ".join(equation_parts) + f" = {np.sum(A[0,:] * x)} = 0 (mod {q})")

    else:
        print("\nSolution is incorrect or the algorithm failed for this specific basis choice.")


if __name__ == "__main__":
    main()
