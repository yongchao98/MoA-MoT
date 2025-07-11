import numpy as np
from numpy.linalg import matrix_power, inv, matrix_rank

def solve():
    """
    Designs a deadbeat observer for the given discrete-time system.
    """
    # Define the system matrices
    A = np.array([
        [-1, 0, 0, 1],
        [1, 0, 0, 2],
        [0, 1, 0, -1],
        [-1, 0, 1, -1]
    ])

    C = np.array([
        [1, 0, 0, 0],
        [1, 0, 0, 1]
    ])

    n = A.shape[0]

    # Step 1: Formulate the dual control problem
    Ad = A.T
    Bd = C.T

    # Step 2: Reduce the multi-input system to a single-input system
    # We choose g = [0, 1]^T. This effectively selects the second column of Bd.
    g = np.array([[0], [1]])
    b_new = Bd @ g

    # Check for controllability of the single-input system
    ControllabilityMatrix = np.hstack([matrix_power(Ad, i) @ b_new for i in range(n)])
    if matrix_rank(ControllabilityMatrix) < n:
        print("The selected single-input system is not controllable. A different 'g' vector should be chosen.")
        return

    # Step 3: Use Ackermann's formula to find the controller gain k_tilde
    # For a deadbeat controller, all desired poles are at 0.
    # The desired characteristic polynomial is p(z) = z^n.
    # So, p(Ad) is simply Ad^n.
    p_Ad = matrix_power(Ad, n)

    # Ackermann's formula: k = [0, ..., 1] * inv(ControllabilityMatrix) * p(A)
    e_n = np.zeros((1, n))
    e_n[0, n - 1] = 1
    
    inv_Co = inv(ControllabilityMatrix)
    k_tilde = e_n @ inv_Co @ p_Ad

    # Step 4: Calculate the observer gain L
    # The controller gain for the dual multi-input system is K = g * k_tilde
    K = g @ k_tilde
    
    # The observer gain is the transpose of K
    L = K.T
    
    # Set small values to zero for cleaner output
    L[np.abs(L) < 1e-10] = 0
    
    print("The designed observer gain matrix L is:")
    print(L)

solve()