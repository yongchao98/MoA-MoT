import numpy as np

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

# Based on the pole placement design for the dual system, we found the
# state-feedback gain K for the dual system (A^T, C^T).
# We chose a solution where the second row of K is zero for simplicity.
K = np.array([
    [-2.0, -1.0, 3.0, -1.0],
    [ 0.0,  0.0, 0.0,  0.0]
])

# The observer gain L is the transpose of K.
L = K.T

print("The designed observer gain matrix L is:")
print(L)

# --- Verification ---
# For a deadbeat observer, all eigenvalues of the error dynamics matrix (A - LC)
# must be at the origin (0).
A_cl = A - L @ C

# Calculate the eigenvalues of the closed-loop system matrix
eigenvalues = np.linalg.eigvals(A_cl)

print("\nTo verify the design, we calculate the eigenvalues of (A - LC).")
print("For a deadbeat observer, all eigenvalues should be zero.")
print("Calculated eigenvalues (rounded to 5 decimal places):")
# We round the result to handle potential floating-point inaccuracies.
print(np.round(eigenvalues, 5))
