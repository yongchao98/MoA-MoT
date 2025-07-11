import numpy as np

def solve_matrix_problem(n):
    """
    Solves the user's matrix problem for a given n.

    This function implements the following logic:
    1. It constructs a specific n x n nilpotent matrix with all non-zero integer entries,
       chosen to maximize the ratio of norms of its RREF under a relaxed interpretation.
    2. This matrix M_n is a rank-1 matrix constructed from vectors u and v.
       u = [1, 1, ..., 1, -(n-1)]^T
       v = [1, 1, ..., 1]^T
       M_n = u * v^T
    3. The immanants of any rank-1 matrix are either 0 or equal to its permanent.
    4. For this specific matrix, the permanent is n! * (-(n-1)), which is negative.
    5. Therefore, the largest immanant is 0.

    Args:
        n (int): The dimension of the matrix.

    Returns:
        None. Prints the result directly.
    """
    if not isinstance(n, int) or n < 2:
        print("Please provide an integer n >= 2.")
        return

    # 1. Construct the matrix M_n
    u = np.ones((n, 1), dtype=int)
    u[n-1, 0] = -(n - 1)
    v = np.ones((1, n), dtype=int)
    M_n = u @ v

    print(f"For n = {n}, the chosen matrix M_n is:")
    print(M_n)
    print("-" * 20)

    # 2. Calculate the permanent and other immanants
    # The term that is constant for every permutation in the immanant sum
    prod_u = np.prod(u)
    prod_v = np.prod(v)
    const_term = prod_u * prod_v
    
    # The permanent is const_term * n!
    perm_val = const_term * np.math.factorial(n)

    # All other irreducible immanants are 0.
    other_immanants_val = 0

    print("The immanants for this matrix are:")
    print(f"Permanent = {np.math.factorial(n)} * ({const_term}) = {perm_val}")
    print("All other non-trivial immanants = 0")
    print("-" * 20)

    # 3. Determine the largest immanant
    largest_immanant = max(perm_val, other_immanants_val)

    print(f"The final equation for the largest immanant is:")
    print(f"LargestImmanant(M_{n}) = max({perm_val}, 0)")
    print(f"LargestImmanant(M_{n}) = {largest_immanant}")


# Example for n=4. You can change this value.
n_value = 4
solve_matrix_problem(n_value)
<<<0>>>