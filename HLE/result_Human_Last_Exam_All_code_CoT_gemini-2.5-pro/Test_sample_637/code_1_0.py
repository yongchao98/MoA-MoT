import numpy as np
import control as ct

def solve_deadbeat_observer():
    """
    This function designs a deadbeat observer for the given discrete-time system
    by calculating the observer gain matrix L.
    """
    # Define the system matrices from the problem description
    A = np.array([[-1, 0, 0, 1],
                  [ 1, 0, 0, 2],
                  [ 0, 1, 0, -1],
                  [-1, 0, 1, -1]])

    C = np.array([[1, 0, 0, 0],
                  [1, 0, 0, 1]])

    # The state-space dimension is n=4.
    n = A.shape[0]
    
    # For a deadbeat observer, all observer poles must be at the origin.
    desired_poles = np.zeros(n)

    # First, we check if the system is observable. Pole placement is only
    # possible for observable systems.
    observability_matrix = ct.obsv(A, C)
    if np.linalg.matrix_rank(observability_matrix) < n:
        print("The system is not observable, so arbitrary pole placement is not possible.")
        return

    # We solve the dual state-feedback controller problem.
    # The dual system is (A.T, C.T). We want to find a gain K such that
    # the eigenvalues of (A.T - C.T * K) are the desired poles.
    # The control.place function finds K for (A - B*K), so we use A.T as the
    # system matrix and C.T as the input matrix.
    try:
        K = ct.place(A.T, C.T, desired_poles)

        # The observer gain L is the transpose of the controller gain K.
        L = K.T

        print("The deadbeat observer gain matrix L is:")
        # Print each row of the matrix for clarity
        for i, row in enumerate(L):
            print(f"  Row {i+1}: {row}")

    except Exception as e:
        print(f"An error occurred during pole placement: {e}")

# Run the function to get the solution
solve_deadbeat_observer()
