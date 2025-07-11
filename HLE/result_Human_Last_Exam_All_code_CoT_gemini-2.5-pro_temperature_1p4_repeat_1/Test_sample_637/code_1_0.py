import numpy as np
import control as ct

# Set print options for a clear matrix output
np.set_printoptions(precision=4, suppress=True)

# Define the system matrices as given in the problem
A = np.array([[-1, 0, 0, 1],
              [1, 0, 0, 2],
              [0, 1, 0, -1],
              [-1, 0, 1, -1]])

C = np.array([[1, 0, 0, 0],
              [1, 0, 0, 1]])

# The number of states determines the number of poles to place
n = A.shape[0]

# For a deadbeat observer, all desired poles are at the origin
poles = np.zeros(n)

# We can design the observer by using pole placement on the dual system.
# The `control.place` function calculates a state-feedback gain K for a system.
# We want to find L such that eigenvalues of (A - LC) are at the origin.
# This is equivalent to finding K = L^T for the dual system (A^T, C^T)
# such that eigenvalues of (A^T - C^T*K) are at the origin.

try:
    # Calculate the state-feedback gain for the dual system (A.T, C.T)
    K_for_dual = ct.place(A.T, C.T, poles)

    # The observer gain L is the transpose of the calculated gain K.
    L = K_for_dual.T

    # Print the final observer gain matrix
    print("The designed observer gain matrix L is:")
    print(L)

except Exception as e:
    print(f"An error occurred during pole placement: {e}")
    print("Please ensure the system (A, C) is observable.")
