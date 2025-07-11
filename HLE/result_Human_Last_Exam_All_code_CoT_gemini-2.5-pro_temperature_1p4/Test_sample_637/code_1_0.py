import numpy as np
from scipy.signal import place_poles

# Define the system matrices from the problem statement
A = np.array([[-1, 0, 0, 1],
              [1, 0, 0, 2],
              [0, 1, 0, -1],
              [-1, 0, 1, -1]])

C = np.array([[1, 0, 0, 0],
              [1, 0, 0, 1]])

# For a deadbeat observer, all observer poles must be placed at the origin.
# The number of poles is equal to the number of states (which is 4).
n_states = A.shape[0]
desired_poles = np.zeros(n_states)

# To find the observer gain L for the poles of (A - LC), we solve the dual
# pole placement problem. We want to find a gain K_dual for the dual system
# defined by A_dual = A.T and B_dual = C.T.
# The poles of the observer will be the poles of (A_dual - B_dual * K_dual).
A_dual = A.T
B_dual = C.T

# Use the `place_poles` function to calculate the dual gain matrix K_dual.
# This function returns a result object containing the gain matrix.
fsf = place_poles(A_dual, B_dual, desired_poles)
K_dual = fsf.gain_matrix

# The desired observer gain L is the transpose of the calculated dual gain.
L = K_dual.T

# Print the final observer gain matrix
print("The observer gain matrix L is:")
print(L)
