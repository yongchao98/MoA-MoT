import numpy as np

def reduce_to_observer_canonical():
    """
    Reduces a discrete-time system to observer canonical form using duality.
    """
    # Define the original system matrices
    A = np.array([[1, 1, 0],
                  [2, 1, 1],
                  [0, 2, 0]])

    C = np.array([[0, 1, 0],
                  [1, 1, 0]])

    n = A.shape[0]

    # --- Step 2: Handle Multi-Output ---
    # Select the first row of C to create a SISO system for the transformation
    # We must ensure the pair (A, c1) is observable.
    c1 = C[0, :]
    observability_matrix_siso = np.vstack([c1 @ np.linalg.matrix_power(A, i) for i in range(n)])
    if np.linalg.matrix_rank(observability_matrix_siso) < n:
        print("The pair (A, c1) is not observable. Please choose a different output row.")
        return

    # --- Step 3: Transform the Dual System ---
    # Form the dual system (Ad, Bd)
    Ad = A.T
    Bd = c1.reshape(-1, 1)

    # Calculate the controllability matrix of the dual system
    Wc = np.hstack([np.linalg.matrix_power(Ad, i) @ Bd for i in range(n)])

    # Find the characteristic polynomial of A: p(s) = s^n + a1*s^(n-1) + ... + an
    # numpy.poly returns [1, a1, a2, ..., an]
    poly_coeffs = np.poly(A)
    a_coeffs = poly_coeffs[1:]

    # Construct the transformation matrix M (W_ccf_inv)
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(n - i):
            if j == 0:
                M[i, j] = a_coeffs[n - 1 - i]
            elif j == 1:
                M[i, j] = a_coeffs[n - 2 - i] if (n - 2 - i) >= 0 else 1
            else:
                if (n - 1 - i - j) >= 0:
                     M[i,j] = a_coeffs[n-1-i-j]
    M = np.fliplr(np.eye(n))
    for i in range(1, n):
        M = np.vstack([a_coeffs[n-i-1:n-1], M[:-1,:]])
    
    # Correct M matrix construction for the form p(s) = s^3+a1s^2+a2s+a3
    # M = [[a2, a1, 1], [a1, 1, 0], [1, 0, 0]]
    M = np.array([
        [a_coeffs[1], a_coeffs[0], 1],
        [a_coeffs[0], 1, 0],
        [1, 0, 0]
    ])
    
    # Calculate P, the transformation matrix for the controller form of the dual
    P = Wc @ M

    # --- Step 4: Find the Observer Transformation T ---
    P_inv = np.linalg.inv(P)
    T = P_inv.T

    # --- Step 5: Calculate New C Matrix ---
    C_new = C @ T

    # Print the results
    print("The new C matrix, C_o, is the result of the matrix multiplication C_o = C * T\n")
    
    np.set_printoptions(precision=4, suppress=True)

    print("Original C matrix:")
    print(C)
    print("\nTransformation matrix T:")
    print(T)
    print("\nResulting new C matrix (in observer canonical form):")
    print(C_new)
    
    # The final answer to be extracted
    print("\n<<<")
    print(C_new)
    print(">>>")


reduce_to_observer_canonical()