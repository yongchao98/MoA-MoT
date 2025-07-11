import numpy as np
import scipy.signal

def solve_deadbeat_observer():
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

    # Get the dimension of the state space
    n = A.shape[0]

    # Step 1: Check for observability
    # Construct the observability matrix O = [C; CA; CA^2; ...; CA^(n-1)]
    O = C
    A_power_i = A
    for _ in range(1, n):
        O = np.vstack([O, C @ A_power_i])
        A_power_i = A_power_i @ A

    # Calculate the rank of the observability matrix
    rank_O = np.linalg.matrix_rank(O)
    
    if rank_O < n:
        print(f"The system is not fully observable (Rank of observability matrix is {rank_O} but n={n}).")
        print("A deadbeat observer that drives the full state error to zero might not exist.")
        # We proceed assuming the unobservable modes are benign (at the origin),
        # which allows for the design. The standard tools should still work if a solution exists.
        
    # Step 2: Define the desired poles for the deadbeat observer
    # All poles should be at the origin for the error to decay to zero in a finite number of steps.
    desired_poles = np.zeros(n)

    # Step 3: Solve the pole placement problem for the dual system (A.T, C.T)
    # The observer gain L is the transpose of the state-feedback gain K for the dual system.
    # The function place_poles finds K such that eigenvalues of (A.T - C.T @ K) are the desired_poles.
    A_T = A.T
    C_T = C.T
    
    try:
        result = scipy.signal.place_poles(A_T, C_T, desired_poles)
        K = result.gain_matrix
    except Exception as e:
        print(f"Error during pole placement: {e}")
        print("This may happen if the system (A.T, C.T) is not controllable, which means (A, C) is not observable.")
        return

    # Step 4: The observer gain L is the transpose of K.
    L = K.T
    
    # Print the final result
    print("The observer gain matrix L is:")
    # Set print options for better readability
    np.set_printoptions(precision=4, suppress=True)
    print(L)

solve_deadbeat_observer()