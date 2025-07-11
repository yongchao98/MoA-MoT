import numpy as np

def solve_and_demonstrate():
    """
    This function demonstrates the solution by constructing a concrete example.
    It calculates the minimal-norm perturbation E for a given A, b, and x,
    and shows that its rank is 1.
    """
    # Set print options for better readability of matrices and vectors
    np.set_printoptions(precision=4, suppress=True)

    # --- Setup of the Example ---
    # We choose A, x, and b such that x is not already a least-squares solution for (A,b).
    # This ensures that the required perturbation E will be non-zero.
    A = np.array([[1.0, 0.0],
                  [0.0, 0.0]])
    x = np.array([[1.0],
                  [0.0]])
    b = np.array([[0.0],
                  [1.0]])

    print("--- Step-by-step Calculation for a Concrete Example ---")
    print("This code will find the minimal Frobenius norm matrix E such that x is a least-squares")
    print("solution to min_z ||(A+E)z - b||_2, and then find its rank.")
    print("\nGiven A:\n", A)
    print("\nGiven x:\n", x)
    print("\nGiven b:\n", b)
    print("-" * 40)

    # Theory recap: The minimal norm E has the form E = c @ x.T
    # where c = - (pinv(A.T * ||x||^2 + r @ x.T)) @ r
    # and r = b - A @ x.

    # 1. Calculate the residual r = b - A @ x
    r = b - A @ x
    print("1. Calculate the residual for the original system: r = b - A @ x")
    print(f"   r = {b.flatten()} - { (A@x).flatten() } = {r.flatten()}")
    print("   r:\n", r)
    print("-" * 40)

    # 2. Calculate the squared norm of x
    x_norm_sq = (x.T @ x).item()
    print(f"2. Calculate the squared norm of x: ||x||^2 = {x_norm_sq:.4f}")
    print("-" * 40)

    # 3. Form the intermediate matrix M = A.T * ||x||^2 + r @ x.T
    M = A.T * x_norm_sq + r @ x.T
    print("3. Form the matrix M = A.T * ||x||^2 + r @ x.T")
    print(f"   M = \n{A.T} * {x_norm_sq:.4f} + \n{r @ x.T}")
    print("   M:\n", M)
    print("-" * 40)

    # 4. Calculate the Moore-Penrose pseudoinverse of M
    M_pinv = np.linalg.pinv(M)
    print("4. Calculate the pseudoinverse of M, M_pinv:")
    print("   M_pinv:\n", M_pinv)
    print("-" * 40)

    # 5. Calculate the vector c = -M_pinv @ r
    c = -M_pinv @ r
    print("5. Calculate the vector c = -M_pinv @ r:")
    print(f"   c = - \n{M_pinv} @ \n{r}")
    print("   c:\n", c)
    print("-" * 40)

    # 6. Calculate E = c @ x.T (the outer product)
    E = c @ x.T
    print("6. Calculate E = c @ x.T:")
    print(f"   E = \n{c} @ \n{x.T}")
    print("   E:\n", E)
    print("-" * 40)

    # 7. Calculate and print the rank of E
    rank_E = np.linalg.matrix_rank(E)
    print(f"7. The rank of the resulting matrix E is {rank_E}.")
    print("-" * 40)

    print("Conclusion:")
    print("This example demonstrates a case where the rank is 1.")
    print("Since the theory shows the rank of E can be at most 1, the greatest possible rank is 1.")

solve_and_demonstrate()
<<<1>>>