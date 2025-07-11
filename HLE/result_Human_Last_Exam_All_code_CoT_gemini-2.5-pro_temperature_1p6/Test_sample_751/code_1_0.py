import numpy as np

# Step 1: Define the system matrices
A = np.array([[1, 1, 0],
              [2, 1, 1],
              [0, 2, 0]])

C = np.array([[0, 1, 0],
              [1, 1, 0]])

# Choose the first row of C to create a single-output system for the transformation
c1 = C[0, :].reshape(1, -1)
print("Original A matrix:\n", A)
print("\nOriginal C matrix:\n", C)
print("\nSelected output vector c1 for transformation:\n", c1)

# Step 2: Form the observability matrix for (A, c1) and check for observability
O1 = np.vstack([c1,
                c1 @ A,
                c1 @ A @ A])

print("\nObservability matrix O1 for (A, c1):\n", O1)

# Check if the system is observable by checking the determinant of O1
det_O1 = np.linalg.det(O1)
print(f"\nDeterminant of O1: {det_O1:.2f}")
if abs(det_O1) < 1e-9:
    print("The system (A, c1) is not observable. Cannot convert to observer canonical form.")
else:
    print("The system (A, c1) is observable.")

    # Step 3: Compute the transformation matrix P
    # The first column p1 is found by solving O1 * p1 = [0, 0, 1]^T
    # This is equivalent to p1 = inv(O1) * [0, 0, 1]^T
    O1_inv = np.linalg.inv(O1)
    e3 = np.array([[0], [0], [1]])
    
    p1 = O1_inv @ e3
    
    # The other columns are found by p_k+1 = A * p_k
    p2 = A @ p1
    p3 = A @ p2
    
    # Construct the transformation matrix P from its columns
    P = np.hstack([p1, p2, p3])
    
    print("\nTransformation matrix P:\n", P)

    # Step 4: Calculate the new C matrix in the observer canonical form basis
    C_new = C @ P
    
    print("\nThe new matrix C in observer canonical form is:")
    print(C_new)

    # For verification, we can also compute the new A matrix (A_o)
    # A_new = np.linalg.inv(P) @ A @ P
    # print("\nNew A matrix (for verification):\n", np.round(A_new, 5))


<<<
[[0.  , 0.  , 1.  ],
 [0.5 , 0.5 , 1.5 ]]
>>>