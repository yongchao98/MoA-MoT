import numpy as np
from scipy import signal

def solve_deadbeat_observer():
    """
    Designs a deadbeat observer for a given discrete-time system.

    The observer gain L is calculated such that all eigenvalues of (A - LC) are zero,
    ensuring the observation error converges to zero in a finite number of steps.
    This is achieved by solving the dual pole placement problem for the system (A^T, C^T).
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

    # The dimension of the state space (n)
    n = A.shape[0]

    # For a deadbeat observer, all eigenvalues of (A - LC) should be at the origin.
    # These are the desired poles for the observer error dynamics.
    desired_poles = np.zeros(n)

    # We solve the dual problem: find a controller gain K for the system (A^T, C^T)
    # such that the poles of (A^T - C^T K) are at the desired locations.
    # The scipy.signal.place_poles function can compute this gain K.
    # The observer gain L is then the transpose of K.
    
    # The dual system matrices are A.T and C.T
    A_dual = A.T
    B_dual = C.T # This is the "B" matrix for the dual control problem

    # Use place_poles to find the state-feedback gain K for the dual system
    # The function returns a full result object, from which we extract the gain matrix.
    placement_result = signal.place_poles(A_dual, B_dual, desired_poles)
    K = placement_result.gain_matrix

    # The observer gain L is the transpose of K
    L = K.T

    # Print the resulting observer gain matrix L
    print("The observer gain matrix L is:")
    print(L)

solve_deadbeat_observer()