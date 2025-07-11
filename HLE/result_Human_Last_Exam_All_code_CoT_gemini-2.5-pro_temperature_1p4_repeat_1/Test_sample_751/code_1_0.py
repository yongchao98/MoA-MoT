import numpy as np

# Step 1: Define the system matrices
A = np.array([
    [1, 1, 0],
    [2, 1, 1],
    [0, 2, 0]
])

C = np.array([
    [0, 1, 0],
    [1, 1, 0]
])

# Step 2: Define the dual system
A_d = A.T
B_d = C.T

print("Original A matrix:\n", A)
print("\nOriginal C matrix:\n", C)
print("\nDual matrix A_d = A^T:\n", A_d)
print("\nDual matrix B_d = C^T:\n", B_d)

# Step 3: Construct the controllability matrix for the dual system
# We use the first input of the dual system (first column of B_d)
b1 = B_d[:, 0].reshape(-1, 1)
b2 = A_d @ b1
b3 = A_d @ b2

Q_c = np.hstack([b1, b2, b3])
print("\nControllability matrix of the dual system (Q_c):\n", Q_c)

# Check if the system is controllable
if np.linalg.det(Q_c) == 0:
    print("\nThe dual system is not controllable with the first input. Cannot proceed.")
else:
    # Step 4: Construct the transformation matrix T
    Q_c_inv = np.linalg.inv(Q_c)
    print("\nInverse of the controllability matrix (Q_c_inv):\n", Q_c_inv)
    
    # The last row of Q_c_inv forms the last row of the transformation matrix T
    t3 = Q_c_inv[-1, :]
    
    # The other rows are found by t_{i-1} = t_i * A_d
    t2 = t3 @ A_d
    t1 = t2 @ A_d
    
    T = np.vstack([t1, t2, t3])
    print("\nTransformation matrix T:\n", T)

    # Step 5: Calculate the new C matrix for the observer canonical form
    # C_o = C * P^-1 = C * T^T
    C_o = C @ T.T
    
    print("\nThe new matrix C in observer canonical form is C_o = C @ T^T:")
    # Print the final matrix C_o, which is the answer
    print(C_o)
