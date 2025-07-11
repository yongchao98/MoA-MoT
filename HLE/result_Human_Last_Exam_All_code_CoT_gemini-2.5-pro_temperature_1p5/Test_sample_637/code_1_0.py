import numpy as np
from scipy.signal import place_poles

# Define the system matrices A and C based on the problem description.
# The system is x(k+1) = Ax(k) + Bu(k), y(k) = Cx(k).
# We are designing an observer, so the B matrix is not needed for the gain calculation.
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

# The state vector has dimension n.
n = A.shape[0]

# For a deadbeat observer, the error must go to zero in the shortest number of steps.
# This requires placing all eigenvalues (poles) of the error dynamics matrix (A - LC) at the origin.
desired_poles = np.zeros(n)

# We can find the observer gain L using pole placement on the dual system.
# The eigenvalues of (A - LC) are the same as the eigenvalues of (A - LC)^T = A^T - C^T * L^T.
# We can use a standard pole placement function to find K_dual = L^T for the system (A^T, C^T).
# The place_poles function finds K such that the eigenvalues of (A - BK) are p.
# For our dual problem, we use A.T for A, C.T for B, and K_dual for K.
try:
    placement_result = place_poles(A.T, C.T, desired_poles)
    K_dual = placement_result.gain

    # The observer gain L is the transpose of the calculated gain K_dual for the dual system.
    L = K_dual.T

    # Print the resulting observer gain matrix L.
    print("The gain matrix L for the deadbeat observer is:")
    print(L)

except Exception as e:
    print("An error occurred during pole placement. This might be because the system is not observable.")
    print(e)
    # Check for observability as a likely cause.
    O_list = [C]
    for i in range(1, n):
        O_list.append(O_list[-1] @ A)
    O = np.vstack(O_list)
    rank = np.linalg.matrix_rank(O)
    print(f"System order n = {n}")
    print(f"Observability matrix rank = {rank}")
    if rank < n:
        print("The system is not observable, so arbitrary pole placement is not possible.")
