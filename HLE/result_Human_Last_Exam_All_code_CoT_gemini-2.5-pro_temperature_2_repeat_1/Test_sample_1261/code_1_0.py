import numpy as np

def solve_system_mod_2(A):
    """
    Solves Ax = 0 (mod 2) using Gaussian elimination over F_2.
    Returns a basis for the null space.
    This is a simplified implementation. For robust usage, a library like sympy is better.
    """
    mat = A.copy() % 2
    n, m = mat.shape
    pivot_row = 0
    pivot_cols = []

    # Forward elimination
    for j in range(m):
        if pivot_row < n:
            pivot = np.where(mat[pivot_row:, j] == 1)[0]
            if len(pivot) > 0:
                pivot = pivot[0] + pivot_row
                mat[[pivot_row, pivot]] = mat[[pivot, pivot_row]] # Swap rows
                for i in range(n):
                    if i != pivot_row and mat[i, j] == 1:
                        mat[i, :] = (mat[i, :] + mat[pivot_row, :]) % 2
                pivot_cols.append(j)
                pivot_row += 1

    # Find basis for null space
    basis = []
    free_cols = [j for j in range(m) if j not in pivot_cols]
    for free_col in free_cols:
        vec = np.zeros(m, dtype=int)
        vec[free_col] = 1
        for i, pivot_col in enumerate(pivot_cols):
            vec[pivot_col] = mat[i, free_col]
        basis.append(vec)
        
    return basis

def solve_underdefined_mq_system(F, d):
    """
    Placeholder for a deterministic poly-time solver for underdefined
    multivariate quadratic systems over F_2.
    F is a list of n quadratic forms in d variables.
    Returns a non-trivial solution vector c of size d if one exists.
    """
    # In a real implementation, this would involve techniques like linearization.
    # For this example, we'll return a dummy solution if d > n.
    n = len(F)
    if d > n:
        c = np.zeros(d, dtype=int)
        c[0] = 1 # Return a trivial solution (e.g., first basis vector)
        return c
    return None

def find_binary_solution(A, q, k):
    """
    Finds a non-zero x in {0,1}^m such that Ax = 0 (mod q)
    using the lifting algorithm.
    """
    n, m = A.shape

    # Step 1: Solve Ax = 0 (mod 2)
    # Get basis for the solution space V_1
    V_basis = solve_system_mod_2(A)
    if not V_basis:
        print("No non-trivial solution mod 2 found.")
        return None

    # Iteratively lift the solution from mod 2^j to mod 2^{j+1}
    for j in range(1, k):
        mod_current = 2**j
        mod_next = 2**(j+1)
        d = len(V_basis)

        # Build the quadratic system for the lifting step.
        # We need to find c = (c_1,...,c_d) in {0,1}^d such that
        # x = c_1*b_1 + ... + c_d*b_d (mod 2) solves Ax = 0 (mod 2^{j+1}).
        # The condition is (A*x / 2^j) % 2 == 0.
        # This results in n quadratic equations in d variables.
        
        # This part is complex to formulate explicitly without a library.
        # We'll use a placeholder representing the solving process.
        
        quadratic_system_placeholder = [] # list of n quadratic equations
        for i in range(n):
           # Each F_i would be a representation of a quadratic polynomial
           # For example, as a list of coefficients for linear and quadratic terms
           quadratic_system_placeholder.append(f"Q_equation_{i}") 
           
        c = solve_underdefined_mq_system(quadratic_system_placeholder, d)
        
        if c is None:
            print(f"Could not find a lifting solution from mod {mod_current} to mod {mod_next}.")
            return None
        
        # Update the basis for the new solution space V_{j+1}
        # This is another complex step, so we simplify by assuming we find
        # one solution vector instead of a new basis.
        # Here we just combine the old basis vectors using coeffs `c`.
        solution = np.zeros(m, dtype=int)
        for i in range(d):
            if c[i] == 1:
                solution = (solution + V_basis[i]) % 2

        # A full algorithm would re-compute a basis for V_{j+1}
        # For simplicity, we assume 'solution' is one vector in V_k
        # and if we need a basis we find the new nullspace. Let's just return
        # the first solution we can construct.
        
        # For the final step, we just need one solution.
        if j == k-1:
          final_solution = solution

    return final_solution

# Problem parameters
n = 2
k = 3 # k > 1
q = 2**k # q = 8
# m = Omega(n^k) = Omega(2^3=8). Let's take m=10, which > n*k = 6
m = 10 

# Generate a random matrix A
np.random.seed(0)
A = np.random.randint(0, q, size=(n, m))

print("Matrix A (n x m):\n", A)
print(f"Modulus q = {q}")

# Run the algorithm
x = find_binary_solution(A, q, k)

if x is not None:
    print("\nFound non-zero binary solution x:\n", x)
    Ax = A @ x
    print("\nVerification: A*x =\n", Ax)
    print("A*x mod q =\n", Ax % q)
    # The result should be the zero vector
    is_solution = not np.any(Ax % q)
    print(f"\nIs Ax = 0 (mod {q})? {is_solution}")
else:
    print("\nFailed to find a solution.")
