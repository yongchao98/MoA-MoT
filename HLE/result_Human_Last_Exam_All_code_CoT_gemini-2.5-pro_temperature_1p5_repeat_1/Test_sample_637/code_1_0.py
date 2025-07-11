import numpy as np
from scipy import signal

def solve_observer_gain():
    """
    Designs a deadbeat observer by placing poles at the origin using the dual system approach.
    """
    # System matrices from the problem description
    A = np.array([[-1, 0, 0, 1], 
                  [1, 0, 0, 2], 
                  [0, 1, 0, -1], 
                  [-1, 0, 1, -1]])
    
    C = np.array([[1, 0, 0, 0], 
                  [1, 0, 0, 1]])

    # The state dimension
    n = A.shape[0]

    # For a deadbeat observer, all desired poles are at the origin.
    desired_poles = np.zeros(n)

    # The dual system for pole placement
    A_dual = A.T
    B_dual = C.T

    # First, verify that the system is observable (which means the dual is controllable)
    observability_matrix = signal.obsv(A, C)
    if np.linalg.matrix_rank(observability_matrix) < n:
        print("System is not observable. A deadbeat observer cannot be designed to control all states.")
        return

    # Use scipy's pole placement function to find the controller gain K for the dual system.
    # The gain K is such that the eigenvalues of (A_dual - B_dual @ K) are the desired_poles.
    place_poles_result = signal.place_poles(A_dual, B_dual, desired_poles, method='YT')
    K = place_poles_result.gain_matrix

    # The observer gain L is the transpose of the controller gain K.
    L = K.T
    
    print("The observer gain matrix L is:")
    # We use a loop to print the matrix nicely, fulfilling the "output each number" instruction.
    # We round to a few decimal places to avoid printing floating-point noise.
    for row in L:
        print(" ".join(f"{num:9.4f}" for num in row))

    # Verification of the result (optional)
    # The eigenvalues of (A - L @ C) should be at the origin.
    # A_obs = A - L @ C
    # computed_poles = np.linalg.eigvals(A_obs)
    # print("\nVerification: Computed eigenvalues of (A - LC):")
    # print(np.round(computed_poles, 4))


if __name__ == '__main__':
    solve_observer_gain()