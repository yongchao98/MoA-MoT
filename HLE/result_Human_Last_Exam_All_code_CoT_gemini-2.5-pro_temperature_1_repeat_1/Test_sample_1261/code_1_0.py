import numpy as np

def find_binary_null_space_basis(A):
    """
    Finds a basis for the null space of a binary matrix A over F_2 using
    Gaussian elimination.
    Args:
      A: A numpy array with entries in {0, 1}.
    Returns:
      A numpy array whose columns form a basis for the null space of A.
      The basis vectors are themselves binary vectors.
    """
    mat = np.copy(A)
    n, m = mat.shape
    pivot_row = 0
    pivot_cols = []
    
    # Forward elimination to get row echelon form
    for j in range(m):  # Iterate through columns
        if pivot_row < n:
            pivot_idx = -1
            for i in range(pivot_row, n):
                if mat[i, j] == 1:
                    pivot_idx = i
                    break
            
            if pivot_idx != -1:
                # Swap rows to bring pivot to the top of the remaining matrix
                mat[[pivot_row, pivot_idx]] = mat[[pivot_idx, pivot_row]]
                pivot_cols.append(j)
                # Eliminate other 1s in the current column
                for i in range(n):
                    if i != pivot_row and mat[i, j] == 1:
                        mat[i, :] = (mat[i, :] + mat[pivot_row, :]) % 2
                pivot_row += 1

    # Identify free variables/columns
    free_cols = [j for j in range(m) if j not in pivot_cols]
    
    # Construct basis vectors from free columns
    basis = []
    for free_col_idx in free_cols:
        vec = np.zeros(m, dtype=int)
        vec[free_col_idx] = 1
        # Back-substitution to find values for pivot variables
        for i, pivot_col_idx in enumerate(pivot_cols):
            if mat[i, free_col_idx] == 1:
                vec[pivot_col_idx] = 1
        basis.append(vec)
        
    if not basis:
        # Trivial null space {0}
        return np.zeros((m, 0), dtype=int)
        
    return np.array(basis).T

def solve_ax_eq_zero_mod_q(A, q):
    """
    Recursively finds a non-zero binary vector x such that Ax = 0 (mod q),
    where q must be a power of 2.
    """
    # Base case of the recursion: q = 1 (i.e., Ax = 0 mod 1)
    # This is always true for any integer vector. We need a non-zero binary vector.
    if q == 1:
        x = np.zeros(A.shape[1], dtype=int)
        if A.shape[1] > 0:
            x[0] = 1
        return x

    # Recursive Step
    # 1. Find a basis for the null space of A mod 2.
    A_mod_2 = A % 2
    B = find_binary_null_space_basis(A_mod_2)
    
    if B.shape[1] == 0:
        # This implies only the trivial solution exists mod 2, which contradicts
        # the problem's premise that m is large enough to guarantee solutions.
        raise RuntimeError("Could not find a non-trivial null space modulo 2.")

    # 2. Form the new matrix for the recursive call.
    # The product A @ B will have all entries divisible by 2.
    A_new = (A @ B) // 2
    
    # 3. Recurse to solve A_new * z = 0 (mod q/2).
    z = solve_ax_eq_zero_mod_q(A_new, q // 2)
    
    # 4. Lift the solution z back to the original space to get x = Bz.
    # Since B and z are binary, the product is computed over integers and then taken mod 2.
    x = (B @ z) % 2
    
    return x

def main():
    """
    Main function to demonstrate the algorithm.
    """
    # Setup parameters based on the problem statement.
    # Let k > 1.
    k = 4
    # Let q = 2^k.
    q = 2**k
    # Let n be a small integer.
    n = 3
    # m must be large enough, m = Omega(n^k) and m > kn.
    # kn = 4 * 3 = 12. Let's choose m > 12.
    m = 15
    # n^k = 3^4 = 81. For a real example m would need to be >> 81, but for
    # demonstration, m > kn is sufficient to ensure the algorithm runs.
    
    print(f"Solving Ax = 0 (mod q) for a non-zero binary vector x.")
    print(f"Parameters: n={n}, m={m}, k={k}, q={q}")
    
    # Generate a random matrix A with entries in Z_q.
    # We use a seed for reproducibility of the example.
    np.random.seed(0)
    A = np.random.randint(0, q, size=(n, m))
    
    print("\nRandomly generated matrix A:")
    print(A)
    
    # Solve the problem
    try:
        x = solve_ax_eq_zero_mod_q(A, q)
        
        # Verify the solution
        result_vector = (A @ x) % q
        
        print("\nFound a non-zero binary solution x:")
        print(x)
        
        print(f"\nVerification: Ax (mod {q}):")
        print(result_vector)
        
        is_correct = np.all(result_vector == 0)
        is_non_zero = np.any(x != 0)
        
        if is_correct and is_non_zero:
            print("\nVerification successful: The solution is correct and non-zero.")
            print("\nDisplaying the equations for Ax = 0 (mod q):")
            for i in range(n):
                equation_str = []
                sum_val = 0
                for j in range(m):
                    if x[j] == 1:
                        term_val = A[i, j]
                        equation_str.append(f"{term_val}")
                        sum_val += term_val
                
                print(f"Row {i}: {' + '.join(equation_str)} = {sum_val}, and {sum_val} mod {q} = {sum_val % q}")
        else:
            print("\nVerification failed.")
            
    except (RuntimeError, ValueError) as e:
        print(f"\nAn error occurred: {e}")

if __name__ == "__main__":
    main()