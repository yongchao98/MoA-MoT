import numpy as np

def find_ker_mod_2(A_in):
    """
    Finds a basis for the null space of matrix A over the field GF(2).
    
    This function requires the 'galois' library. You can install it using:
    pip install galois

    Args:
        A_in: A numpy array (matrix).

    Returns:
        A matrix U whose columns form a basis for the null space of A mod 2.
        Returns a matrix with zero columns if the null space is trivial.
    """
    try:
        import galois
    except ImportError:
        print("Error: The 'galois' library is required for this script.")
        print("Please install it using: pip install galois")
        return None

    GF2 = galois.GF(2)
    A = GF2(A_in)
    
    # The null_space() method returns the basis vectors as rows.
    null_space_basis_rows = A.null_space()
    
    if null_space_basis_rows.shape[0] == 0:
        # Trivial null space {0}. Return a matrix with 0 columns to reflect d=0.
        return np.empty((A_in.shape[1], 0), dtype=int)

    # We want the basis vectors as columns for our matrix U.
    U = null_space_basis_rows.T
    
    # The view converts it back to a standard numpy array.
    return U.view(np.ndarray)

def solve_binary_sis_power_of_2(A, q):
    """
    Solves Ax = 0 (mod q) for a non-zero vector x in {0,1}^m.
    This function implements the recursive lifting algorithm.

    Args:
        A: The input matrix as a numpy array with integer entries.
        q: The modulus, must be a power of 2.

    Returns:
        A non-zero numpy array x with {0,1} entries, or None if no solution is found.
    """
    # Pre-condition check for q.
    if q <= 0 or (q & (q - 1) != 0 and q != 1):
        raise ValueError("Modulus q must be 1 or a power of 2.")

    # Base Case: q=1
    # Ax = 0 (mod 1) is true for any integer vector x.
    # Return a non-zero {0,1} vector, e.g., the first standard basis vector.
    if q == 1:
        m = A.shape[1]
        if m == 0:
            return None # No variables to form a solution
        x = np.zeros(m, dtype=int)
        x[0] = 1
        return x

    # Recursive Step:
    # 1. Find a basis U for the null space of A mod 2.
    U = find_ker_mod_2(A)
    if U is None: # Library not found
        return None
        
    d = U.shape[1] # Dimension of the null space
    if d == 0:
        # The null space is trivial, so no non-zero solution exists in it.
        return None

    # 2. Compute the new matrix A' = (A @ U) / 2 for the recursive call.
    # A @ U results in a matrix where all entries are divisible by 2.
    A_prime_intermediate = A @ U
    A_prime = A_prime_intermediate // 2

    # 3. Recursively call the solver for the smaller-modulus problem.
    z = solve_binary_sis_power_of_2(A_prime, q // 2)

    if z is None:
        return None # Propagate failure up the recursion chain.

    # 4. Construct the solution for the current level.
    # x = U @ z. Note that U and z are {0,1} vectors, but the resulting
    # x is an integer vector. For random A, it's a {0,1} vector with high prob.
    x = U @ z
    
    return x

if __name__ == '__main__':
    # Problem setup
    n = 4
    k = 3
    q = 2**k # q = 8

    # Set m to be large enough, satisfying m = Omega(n^k)
    # For demonstration, let's take a reasonably large m
    m = n * k + 5 # For example, m=17. n^k = 4^3 = 64. A larger m would be more realistic.

    # Generate a random matrix A from Z_q^{n x m}
    # For a truly uniform distribution, use np.random.randint
    np.random.seed(42) # for reproducibility
    A = np.random.randint(0, q, size=(n, m))
    
    print("Finding a non-zero vector x in {0,1}^m such that Ax = 0 (mod q)")
    print("-" * 60)
    print(f"n = {n}, m = {m}, k = {k}, q = {q}")
    print("Matrix A (first 5 columns):\n", A[:,:5])
    print("-" * 60)

    # Find the solution
    solution_x = solve_binary_sis_power_of_2(A.copy(), q)

    if solution_x is not None:
        print("Solution x found:")
        print(solution_x)
        print("\nVerification:")
        # Verify it's a non-zero {0,1} vector
        is_nonzero = np.any(solution_x != 0)
        is_binary = np.all((solution_x == 0) | (solution_x == 1))
        print(f"Is x non-zero? {is_nonzero}")
        print(f"Is x a {{0,1}} vector? {is_binary}")

        # Verify Ax = 0 (mod q)
        if is_nonzero:
            Ax = (A @ solution_x) % q
            is_correct_solution = np.all(Ax == 0)
            print(f"A @ x (mod q) = {Ax}")
            print(f"Is Ax = 0 (mod q)? {is_correct_solution}")
        else:
            print("Algorithm returned the zero vector, which is not a valid solution.")

    else:
        print("Could not find a solution.")
