import numpy as np

# Step 1: Define the original system matrices A and C
A = np.array([[1, 1, 0], 
              [2, 1, 1], 
              [0, 2, 0]])

C = np.array([[0, 1, 0], 
              [1, 1, 0]])

# Step 2: Select the first output vector for the transformation
c1 = C[0, :].reshape(1, -1) # Shape (1, 3)

# Step 3: Define the dual system
A_d = A.T
# B_d is the transpose of c1, making it a column vector
B_d = c1.T  # Shape (3, 1)

# Step 4: Find the transformation P_c for the controller canonical form of the dual system
# 4a. Calculate the controllability matrix of the dual system: W_c = [B_d, A_d*B_d, A_d^2*B_d]
Wc = np.hstack([B_d, A_d @ B_d, A_d @ A_d @ B_d])

# 4b. Find the characteristic polynomial of A (or A_d): det(sI - A) = s^3 + a1*s^2 + a2*s + a3
# The coefficients are returned as [1, a1, a2, a3]
poly_coeffs = np.poly(A)
a1, a2, a3 = poly_coeffs[1], poly_coeffs[2], poly_coeffs[3]

# 4c. Construct the M matrix for transformation
# This specific structure of M corresponds to the standard controller canonical form
M = np.array([[a2, a1, 1],
              [a1, 1,  0],
              [1,  0,  0]])

# 4d. Calculate the transformation matrix P_c
Pc = Wc @ M

# Step 5: Calculate the transformation matrix T_o for the observer canonical form
# T_o is the inverse of the transpose of P_c
To = np.linalg.inv(Pc.T)

# Step 6: Compute the new C matrix for the system in observer canonical form
C_new = C @ To

# Let's verify the transformation. The new state matrix should be in observer canonical form.
# A_new = inv(To) @ A @ To
# The new output vector corresponding to c1 should be c1_new = c1 @ To = [0, 0, 1]

# Output the result
print("Original A matrix:")
print(A)
print("\nOriginal C matrix:")
print(C)
print("\nCharacteristic Polynomial Coefficients (a1, a2, a3):")
print(f"p(s) = s^3 + ({a1:.2f})s^2 + ({a2:.2f})s + ({a3:.2f})")
print("\nTransformation matrix T_o:")
print(To)
print("\nNew A matrix in observer canonical form (A_new = T_o^-1 * A * T_o):")
A_new = np.linalg.inv(To) @ A @ To
print(A_new)
print("\nFinal new C matrix (C_new = C * T_o):")
# We explicitly format the output to show each number in the final matrix equation
print("C_new =")
print(C_new)
