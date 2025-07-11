import numpy as np

# Step 1: Define the original system matrices
A = np.array([
    [1, 1, 0],
    [2, 1, 1],
    [0, 2, 0]
])

C = np.array([
    [0, 1, 0],
    [1, 1, 0]
])

n = A.shape[0]

# Step 2: Construct the dual system
A_d = A.T
B_d = C.T

# Step 3: Create a single-input dual system
# We select the first column of B_d, which corresponds to the first row of C.
b_d = B_d[:, 0:1] # Keep it as a column vector

# Check if the pair (A_d, b_d) is controllable
Wc_check = np.hstack([b_d, A_d @ b_d, A_d @ A_d @ b_d])
if np.linalg.matrix_rank(Wc_check) < n:
    print("The dual system (A^T, c_1^T) is not controllable.")
    print("Cannot transform to controllable canonical form using this output.")
else:
    # Step 4: Find the transformation matrix T
    # 4a: Calculate the characteristic polynomial of A: p(s) = s^n + a_1*s^(n-1) + ... + a_n
    # np.poly(A) returns [1, a_1, a_2, ..., a_n]
    char_poly_coeffs = np.poly(A)
    # The coefficients we need are a_1, a_2, ..., a_n
    coeffs = char_poly_coeffs[1:]
    
    # 4b: Construct the matrix M for the Ackermann formula transformation
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(n - i):
            if i + j + 1 < n:
                M[i, j] = coeffs[n - i - j - 2]
            elif i + j + 1 == n:
                M[i, j] = 1

    # 4c: Construct the controllability matrix for the dual system (A_d, b_d)
    Wc = np.hstack([np.linalg.matrix_power(A_d, i) @ b_d for i in range(n)])
    # For standard controllable canonical form, we need Wc in reverse order of powers
    Wc = np.fliplr(Wc)
    
    # 4d: Calculate the transformation matrix T
    T = Wc @ M
    
    # Step 5: Get the transformation Q for the original system
    T_inv = np.linalg.inv(T)
    Q = T_inv.T
    
    # Step 6: Calculate the new C matrix in observer canonical form
    C_new = C @ Q

    # Output the result
    print("The original matrix A is:")
    print(A)
    print("\nThe original matrix C is:")
    print(C)
    print("\nThe new matrix C in observer canonical form is:")
    # np.round is used to clean up potential floating point inaccuracies
    print(np.round(C_new, decimals=5))

    # We can also verify the new A matrix, A_o, is in observer canonical form
    A_new = np.linalg.inv(Q) @ A @ Q
    # print("\nFor verification, the new matrix A in observer canonical form is:")
    # print(np.round(A_new, decimals=5))

<<<
[[0.  0.  1. ]
 [0.5 0.5 1.5]]
>>>