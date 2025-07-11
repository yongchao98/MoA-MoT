import numpy as np

# Set print options for better readability of the matrix
np.set_printoptions(precision=4, suppress=True)

# Define the system matrices from the problem description
A = np.array([[-1, 0, 0, 1],
              [1, 0, 0, 2],
              [0, 1, 0, -1],
              [-1, 0, 1, -1]])

C = np.array([[1, 0, 0, 0],
              [1, 0, 0, 1]])

n = A.shape[0] # Dimension of the state space

# --- Step 1 & 2: Simplify to SISO and check observability ---
# We hypothesize that the system is observable using only the second output.
# C_siso corresponds to the second row of C.
C_siso = C[1, :].reshape(1, n)

# Construct the observability matrix for the (A, C_siso) pair
# O = [C_siso; C_siso*A; C_siso*A^2; C_siso*A^3]
O_siso = np.vstack([C_siso @ np.linalg.matrix_power(A, i) for i in range(n)])

# Check if the rank is full, which means the system (A, C_siso) is observable.
if np.linalg.matrix_rank(O_siso) < n:
    print("The system is not observable with only the second output. The method fails.")
else:
    # --- Step 3 & 4: Apply Ackermann's formula for a deadbeat observer ---
    # For a deadbeat observer, the desired characteristic polynomial is p(lambda) = lambda^n.
    # The matrix polynomial p(A) is therefore A^n.
    p_A = np.linalg.matrix_power(A, n)
    
    # The formula is l = p(A) * inv(O_siso) * e_n
    # where e_n is a column vector with 1 in the last position.
    O_siso_inv = np.linalg.inv(O_siso)
    e_n = np.zeros((n, 1))
    e_n[-1] = 1
    
    # Calculate the SISO observer gain vector 'l'
    l = p_A @ O_siso_inv @ e_n
    
    # --- Step 5: Construct the full MIMO observer gain matrix L ---
    # Since we only used the second output for feedback, the first column of L is zero.
    L = np.hstack([np.zeros((n, 1)), l])

    print("The observer gain matrix L is:")
    print(L)

    # --- Verification (Optional) ---
    # The eigenvalues of the closed-loop observer matrix (A - LC) should be all zero.
    # A_cl = A - L @ C
    # eigvals_cl = np.linalg.eigvals(A_cl)
    # print("\nFor verification, the eigenvalues of (A - LC) are:")
    # print(eigvals_cl)