import numpy as np
import control as ct

def design_deadbeat_observer():
    """
    Designs a deadbeat observer for a discrete-time system by placing
    the observer poles at the origin using the duality principle.
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

    # State dimension
    n = A.shape[0]

    # Form the dual system for pole placement
    # A_dual = A.T
    # B_dual = C.T
    A_dual = A.transpose()
    B_dual = C.transpose()

    # Define the desired pole locations for the observer.
    # To ensure numerical stability of the 'place' algorithm,
    # we place the poles very close to zero instead of exactly at zero.
    poles = np.array([0.0, 0.0, 1e-9, -1e-9])

    # Use the 'place' command to find the controller gain K for the dual system.
    # This places the eigenvalues of (A_dual - B_dual*K) at the specified locations.
    K = ct.place(A_dual, B_dual, poles)

    # The observer gain L is the transpose of the controller gain K.
    L = K.transpose()

    # Print the resulting observer gain matrix
    # Rounding for cleaner display as the poles were slightly perturbed
    print(np.round(L, 4))

design_deadbeat_observer()