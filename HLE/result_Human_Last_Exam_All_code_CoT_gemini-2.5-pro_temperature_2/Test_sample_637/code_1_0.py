import numpy as np
import control as ct

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

    # For observer pole placement, we use the dual system (A.T, C.T)
    # This is equivalent to finding a feedback gain K for A' - C'K
    A_dual = A.T
    B_dual = C.T

    # We want to place all observer poles at 0 for a deadbeat response
    desired_poles = [0.0, 0.0, 0.0, 0.0]

    # Use a pole placement algorithm to find the gain for the dual system
    # The place function finds K such that eigenvalues of (A - BK) are desired_poles
    # Here, our A is A_dual and B is B_dual
    K_obs = ct.place(A_dual, B_dual, desired_poles)

    # The observer gain L is the transpose of the calculated gain K_obs
    L = K_obs.T

    print("The observer gain matrix L is:")
    print(L)

solve()