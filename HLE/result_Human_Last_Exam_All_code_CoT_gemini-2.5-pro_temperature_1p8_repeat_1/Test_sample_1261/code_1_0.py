import numpy as np

def find_kernel_basis_gf2(matrix):
    """
    Finds a basis for the null space (kernel) of a matrix over GF(2)
    using Gaussian elimination.
    
    Args:
        matrix: A 2D numpy array with 0/1 entries.
        
    Returns:
        A 2D numpy array where columns form a basis for the kernel.
    """
    mat = np.copy(matrix).astype(int)
    num_rows, num_cols = mat.shape
    
    pivot_row = 0
    pivot_cols = []
    
    # Forward elimination to find row echelon form
    for j in range(num_cols):
        if pivot_row < num_rows:
            pivot = -1
            # Find a pivot (a '1') in the current column
            for i in range(pivot_row, num_rows):
                if mat[i, j] == 1:
                    pivot = i
                    break
            
            if pivot != -1:
                # Swap rows to bring the pivot to the current pivot_row
                mat[[pivot_row, pivot]], mat[[pivot, pivot_row]] = mat[[pivot, pivot_row]].copy(), mat[[pivot_row, pivot]].copy()
                
                # Eliminate other 1s in the same column by adding the pivot row
                for i in range(num_rows):
                    if i != pivot_row and mat[i, j] == 1:
                        mat[i, :] = (mat[i, :] + mat[pivot_row, :]) % 2
                pivot_cols.append(j)
                pivot_row += 1

    free_cols = [j for j in range(num_cols) if j not in pivot_cols]
    
    if not free_cols:
        return np.zeros((num_cols, 0), dtype=int)

    # Construct basis vectors from free variables
    basis = []
    for free_col_idx in free_cols:
        basis_vector = np.zeros(num_cols, dtype=int)
        # Set one free variable to 1, others to 0
        basis_vector[free_col_idx] = 1
        
        # Solve for the pivot variables in terms of the free variable
        for i, pivot_col_idx in enumerate(pivot_cols):
            if mat[i, free_col_idx] == 1:
                basis_vector[pivot_col_idx] = 1
        
        basis.append(basis_vector)

    return np.array(basis).T

def solve_binary_mod_power_of_2(A, q):
    """
    Finds a non-zero binary vector x such that Ax = 0 (mod q),
    where q is a power of 2. This implements an iterative lifting algorithm.
    
    Args:
        A: The n x m matrix of integers.
        q: The modulus, a power of 2 (q = 2^k).
        
    Returns:
        A non-zero m x 1 binary vector x satisfying the condition.
    """
    n, m = A.shape
    if not (q > 0 and (q & (q - 1)) == 0):
        raise ValueError("q must be a power of 2.")

    k = int(np.log2(q))

    # Step 1: Find basis for solutions to Ax = 0 (mod 2)
    A_mod_2 = A % 2
    B_current = find_kernel_basis_gf2(A_mod_2)
    
    if B_current.shape[1] == 0:
        print("Error: Could not find non-trivial solution modulo 2.")
        return None

    # Lift the solution from mod 2^j to mod 2^(j+1) for j=1..k-1
    for j in range(1, k):
        mod_current = 2**j
        d_j = B_current.shape[1]

        # For each basis vector b_i, calculate y_i = (A @ b_i) / 2^j
        Y_list = []
        for i in range(d_j):
            b_i = B_current[:, i]
            y_i = (A.astype(np.int64) @ b_i) // mod_current
            Y_list.append(y_i)
        
        Y = np.array(Y_list).T
        Y_mod_2 = Y % 2

        # Find the kernel of Y_mod_2. C is the change-of-basis matrix.
        C = find_kernel_basis_gf2(Y_mod_2)
        
        if C.shape[1] == 0:
            print(f"Error: Lifting failed at step j={j}.")
            return None 

        # Update the basis for the next level. All operations are GF(2).
        B_current = (B_current @ C) % 2

    # Return the first vector from the final basis as the solution
    return B_current[:, 0].reshape(-1, 1)

# ---- Main Execution ----
# Setup problem parameters according to the description
k_param = 3
q_param = 2**k_param
n_param = 4
# We need m > n*k = 12. Let's choose m=13, which satisfies m=O(n^k) and m=poly(n)
m_param = n_param * k_param + 1

print(f"Let k = {k_param}, so q = 2^k = {q_param}")
print(f"Let n = {n_param}")
print(f"Let m = {m_param} (chosen to satisfy m > nk)")

# Generate a random input matrix A
np.random.seed(42)
A_matrix = np.random.randint(0, q_param, size=(n_param, m_param))
print(f"\nGenerated a random {n_param}x{m_param} matrix A.")

# Solve the problem
x_solution = solve_binary_mod_power_of_2(A_matrix, q_param)

# Print and verify the solution
if x_solution is not None:
    print("\nFound a non-zero binary solution x.")
    
    print("\nThe equation Ax = 0 (mod q) written out with the found solution:")
    for i in range(n_param): # For each row/equation
        line_parts = []
        row_sum = 0
        for j in range(m_param): # For each variable
            term = A_matrix[i, j] * x_solution[j, 0]
            # Show each term in the sum
            line_parts.append(f"({A_matrix[i, j]}*{x_solution[j, 0]})")
            row_sum += term
        
        equation = " + ".join(line_parts)
        print(f"\nEquation for row {i}:")
        print(f"  {equation}")
        print(f"  = {row_sum}")
        print(f"  And {row_sum} mod {q_param} = {row_sum % q_param}")
    
    # Final check
    result_vector = (A_matrix @ x_solution) % q_param
    if np.all(result_vector == 0):
        print("\nVerification successful: Ax is indeed the zero vector mod q.")
else:
    print("\nNo solution was found (this is unexpected given the problem parameters).")
