import numpy as np
from sympy import Matrix

def find_nullspace_mod2(mat):
    """
    Finds a basis for the nullspace of a matrix over the finite field Z_2.
    
    Args:
        mat (np.ndarray): The input matrix with integer entries.
        
    Returns:
        np.ndarray: A matrix whose columns form a basis for the nullspace of mat (mod 2).
                    Returns an empty matrix if the nullspace is trivial or dimensions are zero.
    """
    # If there are no columns (variables), the nullspace basis is empty.
    if mat.shape[1] == 0:
        return np.empty((0, 0), dtype=int)
        
    # Convert numpy array to a sympy Matrix to use its methods
    # We want to find x such that mat * x = 0 (mod 2)
    sympy_mat = Matrix(mat % 2)
    
    # .nullspace() returns a list of vectors that form the basis for the nullspace
    nullspace_basis_vectors = sympy_mat.nullspace()
    
    # If the nullspace is trivial (contains only the zero vector), return an empty basis
    if not nullspace_basis_vectors:
        return np.empty((mat.shape[1], 0), dtype=int)
        
    # Convert the list of sympy vectors to a single numpy matrix
    # Each vector in the list is a column of the basis matrix
    # First, convert each sympy Vector to a list of its elements
    # Then create a numpy array and transpose it to get columns as basis vectors
    basis_matrix = np.array([list(v) for v in nullspace_basis_vectors], dtype=int).T
    
    return basis_matrix

def solve_binary_kernel(A, n, m, k, q):
    """
    Finds a non-zero binary vector x such that Ax = 0 (mod q), where q=2^k.
    Implements the iterative lifting algorithm.
    """
    print("Starting the lifting algorithm...")
    print("-" * 20)

    # Step 1: Find solutions for Ax = 0 (mod 2)
    # The columns of X form a basis for the solution space.
    print("Step 1: Solving Ax = 0 (mod 2)")
    X = find_nullspace_mod2(A)
    print(f"Dimension of solution space mod 2: {X.shape[1]}")

    if X.shape[1] == 0:
        print("No non-trivial binary solution mod 2 found. This is unexpected given the problem constraints.")
        return None

    # Iteratively lift the solution from mod 2^j to mod 2^(j+1)
    for j in range(2, k + 1):
        print("-" * 20)
        print(f"Step {j}: Lifting solution from mod {2**(j-1)} to mod {2**j}")
        
        # We need to solve Mc = 0 (mod 2) where A*X = 2**(j-1) * M
        # Using integer arithmetic for this part.
        T = np.dot(A, X)

        # All entries of T must be divisible by 2**(j-1) by construction.
        if not np.all((T % (2**(j-1))) == 0):
             print(f"Error: Lifting failed at step {j}. Not all elements are divisible by {2**(j-1)}.")
             return None

        M = T // (2**(j-1))
        
        # Find the basis for the nullspace of M (mod 2)
        B = find_nullspace_mod2(M)
        
        print(f"Dimension of new solution space mod {2**j}: {B.shape[1]}")
        
        if B.shape[1] == 0:
            print(f"Lifting failed: no non-trivial solution for the kernel of M at step {j}.")
            return None
            
        # Update the solution basis by right-multiplying with B
        # The new solutions are linear combinations of old solutions
        X = np.dot(X, B) % 2

    # After the loop, columns of X are solutions to Ax = 0 (mod q)
    if X.shape[1] > 0:
        # Return the first basis vector as the solution
        return X[:, 0]
    else:
        return None

if __name__ == '__main__':
    # Problem Parameters
    # Let n, k be small integers > 1
    n = 3
    k = 4
    # m = Omega(n^k) = Omega(3^4) = Omega(81). 
    # Also need m > k*n = 4*3 = 12 to guarantee a solution with this algorithm.
    m = 20 # Using a smaller m for demonstration, it still satisfies m > kn
    q = 2**k
    
    print(f"Parameters: n={n}, k={k}, q=2^{k}={q}, m={m}")
    if m <= k*n:
        print(f"Warning: m={m} is not greater than k*n={k*n}. A solution may not be found by this algorithm.")

    # Generate a random n x m matrix A with entries in Z_q
    A = np.random.randint(0, q, size=(n, m))
    
    # Find the non-zero binary vector x
    x = solve_binary_kernel(A, n, m, k, q)
    
    print("-" * 20)
    
    if x is not None:
        print("Successfully found a non-zero binary solution vector x.")
        print("x =", x)
        
        # Verification
        print("\nVerifying the solution Ax = 0 (mod q)...")
        result_vector = np.dot(A, x) % q
        
        # Output the numbers in the final equation
        for i in range(n):
            # To keep it readable, we show a snippet of the equation
            equation_parts = [f"{A[i, j]}*{x[j]}" for j in range(min(m, 8))]
            ellipsis = " + ..." if m > 8 else ""
            equation_str = " + ".join(equation_parts) + ellipsis
            print(f"Row {i}: ( {equation_str} ) mod {q} = {result_vector[i]}")

        print("\nFinal result vector Ax (mod q):")
        print(result_vector)
        
        if np.all(result_vector == 0):
            print("\nVerification successful: Ax is the zero vector mod q.")
        else:
            print("\nVerification failed: Ax is NOT the zero vector mod q.")
    else:
        print("Could not find a solution.")
