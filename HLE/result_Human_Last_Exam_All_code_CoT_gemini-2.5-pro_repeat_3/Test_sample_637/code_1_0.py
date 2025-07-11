import numpy as np
from scipy.signal import place_poles

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

# The state dimension is n
n = A.shape[0]

# For a deadbeat observer, all eigenvalues of the error dynamics
# matrix (A - LC) must be placed at the origin.
desired_poles = np.zeros(n)

# We solve the dual state-feedback pole placement problem.
# We want to find a gain K for the system (A.T, C.T) such that
# the eigenvalues of (A.T - C.T @ K) are the desired poles.
# The observer gain L is the transpose of this K.
# We use scipy's place_poles function for this.
# result = place_poles(A_sys, B_sys, poles)
# Here, A_sys = A.T and B_sys = C.T
# The gain K will be L.T
try:
    result = place_poles(A.T, C.T, desired_poles)
    K = result.gain_matrix
    
    # The observer gain L is the transpose of the feedback gain K from the dual system.
    L = K.T

    # --- Final Answer Output ---
    # The problem asks to provide the answer as a gain matrix.
    print("The observer gain matrix L is:")
    print(L)

    # --- Verification (Optional) ---
    # To verify, we can compute the eigenvalues of (A - L @ C)
    # and check if they are all close to zero.
    A_obs = A - L @ C
    eigenvalues = np.linalg.eigvals(A_obs)
    # print("\nVerification: Eigenvalues of (A - LC):")
    # print(eigenvalues)

except ValueError as e:
    print(f"Pole placement failed. The system might not be observable. Error: {e}")
