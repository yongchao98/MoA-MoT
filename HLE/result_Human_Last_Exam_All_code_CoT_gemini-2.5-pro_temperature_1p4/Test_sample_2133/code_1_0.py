import numpy as np
from scipy.special import factorial
import math

def solve_matrix_problem(n):
    """
    Solves the user's problem for a given matrix size n.
    
    This function defines a specific 2-nilpotent matrix Mn with non-zero integer
    entries, calculates its Popov normal form (RREF), the ratio of its
    logarithmic mu_infinity norm to its Frobenius norm, and its largest immanant.
    
    Args:
        n (int): The dimension of the square matrix.
    """
    if n < 2:
        print("The matrix dimension n must be at least 2.")
        return

    # Step 1 & 2: Define the matrix M_n
    # We choose a 2-nilpotent matrix with all non-zero integer entries.
    # All rows of the matrix are identical, defined by the vector r.
    r = np.ones(n, dtype=int)
    r[-1] = 1 - n
    
    # Check if any entry is zero, which happens if 1-n=0, i.e., n=1
    # but we already handled n<2. For n>=2, 1-n is a non-zero integer.
    
    M_n = np.tile(r, (n, 1))

    print(f"For n = {n}, the selected matrix M_n is:")
    print(M_n)
    print("-" * 30)

    # Step 3: Compute the Popov normal form (RREF) and the ratio
    # Since all rows are identical, the RREF will have the first row as r
    # and all other rows as zero.
    P_n = np.zeros((n, n), dtype=int)
    P_n[0, :] = r

    print("Its Popov normal form (RREF) P(M_n) is:")
    print(P_n)
    print("-" * 30)
    
    # Calculate the matrix measure mu_infinity(P_n)
    # mu_inf(A) = max_i (A_ii + sum_{j!=i} |A_ij|)
    # For P_n, only the first row is non-zero.
    # For i=0 (first row): p_00 + sum_{j!=0} |p_0j| = r[0] + sum_{j=1}^{n-1} |r[j]|
    # r[0] is 1. The sum is (n-2)*|1| + |1-n| = n-2 + n-1 = 2n-3.
    # So mu_inf for the first row is 1 + 2n-3 = 2n-2.
    # For all other rows, mu_inf is 0.
    # Thus, mu_inf(P_n) is 2n-2.
    mu_inf_norm = 2 * n - 2
    
    # Calculate the Frobenius norm of P_n
    # ||P_n||_F = sqrt(sum of squares of elements)
    # The only non-zero elements are in the first row, which is r.
    frobenius_norm = np.linalg.norm(r)

    # Calculate the ratio
    ratio = mu_inf_norm / frobenius_norm

    print("Logarithmic mu-infinity norm of P(M_n):", mu_inf_norm)
    print("Frobenius norm of P(M_n):", frobenius_norm)
    print(f"Ratio = {mu_inf_norm} / {frobenius_norm} = {ratio}")
    print("-" * 30)
    
    # Step 4: Calculate the largest immanant of M_n
    # For a matrix with identical rows, all immanants corresponding to non-trivial
    # characters of S_n are zero. The only non-zero immanant is the permanent,
    # which corresponds to the trivial character.
    # per(A) = n! * product of the elements of any row, for a matrix A with n identical rows.
    
    # product of elements in row r
    row_product = np.prod(r)
    
    # permanent is n! * row_product
    # math.factorial handles large integers starting Python 3.9
    try:
        # Use math.factorial for potentially large n
        perm = math.factorial(n) * row_product
    except (ValueError, OverflowError):
        # Fallback for very large n if math.factorial fails
        print("n is too large to compute the factorial directly.")
        perm = "too large to compute"
        
    largest_immanant = perm
    
    print(f"The largest immanant of M_n is the permanent.")
    print(f"per(M_n) = {n}! * ({' * '.join(map(str, r))})")
    print(f"per(M_n) = {n}! * {row_product}")
    print(f"Largest immanant value: {largest_immanant}")

# Example for n=3
solve_matrix_problem(3)

# <<<48>>>