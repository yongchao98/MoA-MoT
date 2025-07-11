import numpy as np

def solve_observer_canonical_form():
    """
    Reduces a discrete-time system to observer canonical form using duality.
    """
    # Step 1: Define the original system matrices
    A = np.array([
        [1, 1, 0],
        [2, 1, 1],
        [0, 2, 0]
    ])
    C = np.array([
        [0, 1, 0],
        [1, 1, 0]
    ])

    print("Original A matrix:\n", A)
    print("\nOriginal C matrix:\n", C)

    # Step 2: Define the dual system
    F = A.T
    G = C.T

    print("\nDual system matrix F = A^T:\n", F)
    print("\nDual system matrix G = C^T:\n", G)

    # Step 3: Find the transformation to controller canonical form for the dual system
    # We use the first column of G to build the controllability matrix
    g1 = G[:, 0:1] # Keep it as a column vector

    # Controllability matrix for (F, g1)
    Wc = np.hstack([g1, F @ g1, F @ F @ g1])
    print("\nControllability matrix for the dual system (Wc):\n", Wc)

    # Check if the system is controllable
    if np.linalg.matrix_rank(Wc) < A.shape[0]:
        print("\nDual system is not controllable with the first column of C^T. Cannot proceed.")
        return

    # Characteristic polynomial of F: s^n + a1*s^(n-1) + ... + an
    # np.poly returns coefficients [1, a1, a2, a3]
    coeffs = np.poly(F)
    a1 = coeffs[1]
    a2 = coeffs[2]
    
    # Transformation matrix based on coefficients (Ogata's method)
    Wc_bar = np.array([
        [a2, a1, 1],
        [a1, 1,  0],
        [1,  0,  0]
    ])
    print("\nCoefficient matrix (Wc_bar):\n", Wc_bar)

    # The transformation matrix T for z = T*z_c is T = Wc * Wc_bar
    T_ctrl = Wc @ Wc_bar
    
    # The transformation P_c for z_c = P_c*z is the inverse of T_ctrl
    P_c = np.linalg.inv(T_ctrl)
    print("\nTransformation matrix for dual system controller form (P_c):\n", P_c)

    # Step 4: Use duality to find the new C matrix for the original system
    # The inverse of the observer transformation matrix is Po_inv = P_c.T
    Po_inv = P_c.T
    
    # The new output matrix is Co = C * Po_inv
    C_new = C @ Po_inv
    
    print("\n--- Final Result ---")
    print("The new C matrix in observer canonical form is:")
    print(C_new)

solve_observer_canonical_form()