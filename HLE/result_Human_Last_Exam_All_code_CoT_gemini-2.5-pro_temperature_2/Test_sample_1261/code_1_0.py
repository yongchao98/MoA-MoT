import numpy as np
from sympy import Matrix

def solve_homogeneous_sis(A, q):
    """
    Finds a non-zero vector x in {0,1}^m such that Ax = 0 (mod q).
    
    This function implements a lifting algorithm.
    q is assumed to be a power of 2, q = 2^k.
    A is an n x m matrix.
    The problem constraints on n, m, k are assumed to hold.
    """
    n, m = A.shape
    if not q > 1 or (q & (q-1) != 0):
        raise ValueError("q must be a power of 2 greater than 1.")
        
    k = q.bit_length() - 1
    
    # Step 1: Solve Ax = 0 (mod 2) to get the initial basis.
    # We find the nullspace of A (mod 2) over the field F_2.
    A_mod2 = A % 2
    
    # Use sympy's nullspace over rationals, then convert to a valid F_2 integer basis.
    # This works because for a {0,1} matrix, a rational basis can be cleared of 
    # denominators to yield an integer basis, which is then taken modulo 2.
    sympy_basis = Matrix(A_mod2).nullspace()
    if not sympy_basis:
        print("No non-trivial solution for mod 2. This shouldn't happen with the problem's constraints.")
        return None
    
    # Convert list of sympy Matrices to a numpy array, our basis matrix.
    current_basis_vectors = []
    for vec in sympy_basis:
        # Get lcm of denominators to get an integer vector
        lcm_den = sympy.lcm([f.q for f in vec])
        int_vec = np.array((vec * lcm_den).tolist(), dtype=int).flatten()
        current_basis_vectors.append(int_vec % 2)
    
    # Basis is an m x d_1 matrix, where d_1 is dimension of nullspace.
    current_basis = np.array(current_basis_vectors).T

    # Steps 2 to k: Lift the solution.
    for j in range(1, k):
        p = 2**j # Current power of 2
        
        # d_j is the number of vectors in the current basis
        d_j = current_basis.shape[1]
        
        if d_j == 0:
            print(f"Basis became empty at lifting step j={j}. Cannot find a solution.")
            return None

        # Build the matrix V for the next linear system
        V_cols = []
        for i in range(d_j):
            b_i = current_basis[:, i]
            # Ax is divisible by p=2^j. We want the next bit.
            v_i_exact = A @ b_i
            v_i_lifted = (v_i_exact // p) % 2
            V_cols.append(v_i_lifted)
        
        # V is an n x d_j matrix
        V_matrix = np.array(V_cols).T

        # Find the null space of Vc = 0 (mod 2) to get combination coefficients
        C_sympy_basis = Matrix(V_matrix).nullspace()
        
        if not C_sympy_basis:
            print(f"Lifting failed at step j={j}. No non-trivial combinations found.")
            return None
            
        C_basis_vectors = []
        for vec in C_sympy_basis:
            lcm_den = sympy.lcm([f.q for f in vec])
            int_vec = np.array((vec * lcm_den).tolist(), dtype=int).flatten()
            C_basis_vectors.append(int_vec % 2)

        C_matrix = np.array(C_basis_vectors).T

        # Update the basis for the next level mod 2^{j+1}
        current_basis = (current_basis @ C_matrix) % 2
    
    if current_basis.shape[1] == 0:
        print("Final basis is empty, which implies no non-trivial solution was found.")
        return None
        
    # Any column of the final basis is a solution.
    x = current_basis[:, 0]
    
    # Ensure the returned solution is non-zero
    if np.all(x == 0) and current_basis.shape[1] > 1:
        # If the first basis vector is zero, try the next one
        for i in range(1, current_basis.shape[1]):
            if not np.all(current_basis[:, i] == 0):
                x = current_basis[:, i]
                break
    
    if np.all(x == 0):
        print("Algorithm resulted in a zero vector.")
        return None
        
    return x

def main():
    # Setup problem parameters according to the problem statement
    # k > 1, n >= 1, m=Omega(n^k) and m=poly(n)
    n = 2
    k = 3
    q = 2**k  # q = 8
    
    # We choose m > n*k = 6 for the algorithm to be guaranteed to work
    m = 8
    
    # A is sampled uniformly from Z_q^{n x m}
    np.random.seed(42) # For reproducible results
    A = np.random.randint(0, q, size=(n, m))
    
    print(f"Finding a non-zero binary vector x such that Ax = 0 (mod {q})")
    print("-" * 50)
    print(f"Parameters: n={n}, m={m}, k={k}, q={q}")
    print("\nRandom Matrix A:")
    print(A)

    x = solve_homogeneous_sis(A, q)
    
    if x is not None:
        print("\nFound a non-zero binary solution vector x:")
        print(x)
        
        result_vector = A @ x
        
        print("\nVerification of the equation Ax = 0 (mod q):")
        # Print each number in each equation
        for i in range(n):
            eq_str = " + ".join([f"{A[i, j]}*{x[j]}" for j in range(m)])
            print(f"  Row {i+1}: {eq_str} = {result_vector[i]}")
        
        is_correct = np.all((result_vector % q) == 0)
        print(f"\nResult vector mod {q}: {result_vector % q}. Verification successful: {is_correct}")

    else:
        print("\nFailed to find a solution.")

if __name__ == '__main__':
    main()