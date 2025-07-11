import numpy as np

# Set print options for better readability
np.set_printoptions(precision=4, suppress=True)

# Step 1: Define the original system matrices
A = np.array([[1, 1, 0],
              [2, 1, 1],
              [0, 2, 0]])

C = np.array([[0, 1, 0],
              [1, 1, 0]])

# Step 2: Form the dual system
# The dual system is x_d(k+1) = A_d*x_d(k) + B_d*u_d(k)
# where A_d = A^T and B_d = C^T
A_d = A.T
B_d = C.T

# Step 3: Find the transformation matrix T_c for controller canonical form
# First, we select n=3 linearly independent columns from the controllability matrix
# Co = [B_d, A_d*B_d, ...].
# The columns of B_d are b1 and b2.
b1 = B_d[:, 0]
b2 = B_d[:, 1]
Ad_b1 = A_d @ b1

# The standard procedure for selecting linearly independent vectors is to check
# b1, b2, A_d*b1, A_d*b2, ... in order.
# The set {b1, b2, A_d*b1} is linearly independent for this system.
# We form the matrix M from these vectors.
M = np.column_stack([b1, b2, Ad_b1])

# This selection implies controllability indices mu_1=2 (from the b1 chain)
# and mu_2=1 (from the b2 chain). The sum is 3, the order of the system.

# The transformation matrix T_c is constructed using rows from the inverse of M.
if np.linalg.matrix_rank(M) == A.shape[0]:
    M_inv = np.linalg.inv(M)

    # The transformation T_c depends on rows from M_inv corresponding to the end of each
    # vector chain used to build M. For this problem, these are the 2nd and 3rd rows
    # of M_inv (using 1-based indexing).
    # With mu_1=2 and mu_2=1, the construction rule specifies using the 2nd (sigma_1)
    # and 3rd (sigma_2) rows from the final T_c matrix definition.
    # The generating vectors from M_inv are q_2 (for mu_1 chain) and q_3 (for mu_2 chain).
    # Note: Confusingly, these are the 2nd and 3rd rows of M_inv, which correspond
    # to the last vectors in the construction of the canonical form basis.
    
    q_row_2 = M_inv[1, :] # Corresponds to b2, last in its chain
    q_row_3 = M_inv[2, :] # Corresponds to Ad*b1, last in its chain

    # Construct T_c based on standard multi-input canonical form procedure.
    # Block for mu_1=2 uses q_row_3: t_rows = [q_row_3 @ A_d, q_row_3] -> Incorrect interpretation.
    # Correct (Chen) procedure using M = [b1,b2,Adb1], mu1=2, mu2=1 leads to using
    # M_inv rows corresponding to column 2 (b2) and col 3 (Adb1)
    # The block structure requires using q_sigma_i. sigma_1=mu_1=2, sigma_2=mu_1+mu_2=3.
    # This refers to rows of T, not M_inv directly.
    # The generating rows from M_inv are q_2 and q_3.
    q_gen1 = M_inv[1,:] # q2
    q_gen2 = M_inv[2,:] # q3

    # The rows of Tc are constructed as:
    # {q_gen1 @ A_d, q_gen1} for mu_1=2 block
    # {q_gen2} for mu_2=1 block
    
    # Let's use the consistent calculation from the thinking process
    # M = [b1, b2, Ad*b1] -> inv(M)
    # mu_1=2, mu_2=1 -> sigma_1=2, sigma_2=3
    # Use q2 and q3 from M_inv
    # T_c rows: [q2*Ad, q2, q3]
    q2 = M_inv[1,:] # second row
    q3 = M_inv[2,:] # third row

    t_row1 = q2 @ A_d
    t_row2 = q2
    t_row3 = q3

    T_c = np.vstack([t_row1, t_row2, t_row3])

    # Step 4: Calculate the B matrix in controller canonical form (B_c)
    B_c = T_c @ B_d

    # Step 5: The new C matrix (C_o) is the transpose of B_c
    C_o = B_c.T

    # Print the final result in the specified format
    print("The new matrix C in observer canonical form is C_o = B_c^T.")
    print("Final equation:")
    print(f"C_o = [[{C_o[0, 0]:.4f}, {C_o[0, 1]:.4f}, {C_o[0, 2]:.4f}],")
    print(f"       [{C_o[1, 0]:.4f}, {C_o[1, 1]:.4f}, {C_o[1, 2]:.4f}]]")

else:
    print("The dual system is not controllable, so the original system is not observable.")
    print("Cannot reduce to observer canonical form.")
