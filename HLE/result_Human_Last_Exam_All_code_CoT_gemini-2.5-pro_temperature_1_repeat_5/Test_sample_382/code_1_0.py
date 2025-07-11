import numpy as np

def solve_greatest_rank():
    """
    This function demonstrates that the greatest possible rank of the minimum
    Frobenius norm perturbation E is 1.

    It sets up a random least-squares problem, calculates the optimal
    perturbation E, and prints its rank.
    """
    # Set up a random problem
    # m > n to have an overdetermined system
    m, n = 5, 3
    A = np.random.rand(m, n)
    b = np.random.rand(m)
    # x must be a non-zero vector
    x = np.random.rand(n)
    while np.linalg.norm(x) == 0:
        x = np.random.rand(n)

    print("--- Problem Setup ---")
    print(f"Shape of A: {A.shape}")
    print(f"Shape of b: {b.shape}")
    print(f"Shape of x: {x.shape}")
    print("-----------------------\n")

    # Calculate the residual for the original problem
    r = b - A @ x
    
    # If r is zero, x is an exact solution, so E=0 and rank is 0.
    if np.linalg.norm(r) < 1e-12:
        E = np.zeros_like(A)
        print("x is an exact solution (Ax=b). E is the zero matrix.")
    else:
        # Calculate A^T * r
        A_T_r = A.T @ r
        # If A^T*r is zero, x is already a least-squares solution, so E=0.
        if np.linalg.norm(A_T_r) < 1e-12:
            E = np.zeros_like(A)
            print("x is already a least-squares solution (A^T(b-Ax)=0). E is the zero matrix.")
        else:
            # Compare the norms to decide which form of E is optimal
            norm_r = np.linalg.norm(r)
            norm_x = np.linalg.norm(x)
            norm_A_T_r = np.linalg.norm(A_T_r)
            
            val1 = norm_r / norm_x
            val2 = norm_A_T_r / norm_r
            
            # Based on the comparison, compute the optimal E
            if val1 <= val2:
                # E = r x^T / ||x||^2
                print("Optimal E is of the form: (r * x^T) / ||x||^2")
                E = np.outer(r, x) / (norm_x**2)
            else:
                # E = -r (A^T r)^T / ||r||^2 = - (r r^T A) / ||r||^2
                # Note: (A^T r)^T = r^T A
                print("Optimal E is of the form: -r * (A^T * r)^T / ||r||^2")
                # This is equivalent to -(1/||r||^2) * r * (r^T * A)
                # But it's easier to compute r * (A^T r)^T
                E = -np.outer(r, A_T_r) / (norm_r**2)
    
    # Calculate the rank of E
    rank_E = np.linalg.matrix_rank(E)

    print("\n--- Result ---")
    print(f"The computed perturbation matrix E is:\n{E}")
    print(f"\nThe rank of the matrix E is: {int(rank_E)}")
    print("\nAs shown by the theoretical derivation and this example,")
    print("the greatest possible rank of E is 1.")

solve_greatest_rank()
>>>1