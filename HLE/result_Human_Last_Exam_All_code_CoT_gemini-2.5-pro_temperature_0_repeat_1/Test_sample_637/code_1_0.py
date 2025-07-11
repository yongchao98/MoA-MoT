import numpy as np
import control as ct

def solve_observer_gain():
    """
    Designs a deadbeat observer for the given discrete-time system.
    """
    # Define the system matrices from the problem description
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

    # The dimension of the state vector
    n = A.shape[0]

    # For a deadbeat observer, all poles must be placed at the origin.
    desired_poles = np.zeros(n)

    # We solve the pole placement problem for the dual system.
    # The state matrix of the dual system is A.T (A transpose).
    # The input matrix of the dual system is C.T (C transpose).
    # The control.place function finds a gain K such that the eigenvalues of
    # (A.T - C.T @ K) are placed at the desired locations.
    try:
        K = ct.place(A.T, C.T, desired_poles)
    except Exception as e:
        print(f"An error occurred during pole placement: {e}")
        print("This may happen if the system (A, C) is not observable.")
        return

    # The observer gain L is the transpose of the calculated controller gain K.
    L = K.T

    # Print the final observer gain matrix L
    print("The observer gain matrix L is:")
    print(L)

    # Verification (optional): Check if the eigenvalues of (A - LC) are at 0
    # A_closed_loop = A - L @ C
    # eigenvalues = np.linalg.eigvals(A_closed_loop)
    # print("\nVerification: Eigenvalues of (A - LC):")
    # print(np.round(eigenvalues, 5))


if __name__ == '__main__':
    solve_observer_gain()
