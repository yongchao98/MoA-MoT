import numpy as np

# Given system matrices
A = np.array([
    [1, 1, 0],
    [2, 1, 1],
    [0, 2, 0]
])

C = np.array([
    [0, 1, 0],
    [1, 1, 0]
])

# Step 1: Define the dual system
A_d = A.T
C_d = C.T

# Step 2: Select the first input vector of the dual system
# This corresponds to the first row of the original C matrix
b_d = C_d[:, 0].reshape(-1, 1)

# Determine the dimension of the state space
n = A.shape[0]

# Compute the controllability matrix for the dual system (A_d, b_d)
Wc = b_d
for i in range(1, n):
    Wc = np.hstack((Wc, np.linalg.matrix_power(A_d, i) @ b_d))

# Check if the system is controllable
if np.linalg.matrix_rank(Wc) != n:
    print("The dual system (A^T, c1^T) is not controllable.")
    print("Cannot convert to controller canonical form.")
else:
    # Step 3: Find the transformation matrix T for the dual system
    
    # Get coefficients of the characteristic polynomial of A_d
    # p(lambda) = lambda^n + coeffs[1]*lambda^(n-1) + ... + coeffs[n]
    coeffs = np.poly(A_d)
    
    # Construct the M matrix for the transformation T_inv = Wc * M
    M = np.eye(n)
    for i in range(n):
      for j in range(n - i):
        if i == 0:
            M[j, n-1-j] = 1
        else:
            M[j, n-1-i-j] = coeffs[i]

    # Calculate T_inv
    T_inv = Wc @ M

    # Calculate T by inverting T_inv
    T = np.linalg.inv(T_inv)
    
    # Step 4 & 5: Calculate the new C matrix for the original system
    # The inverse transformation for the original system is P_inv = T^T
    P_inv = T.T
    
    # The new output matrix is C_new = C @ P_inv
    C_new = C @ P_inv
    
    print("The original A matrix is:")
    print(A)
    print("\nThe original C matrix is:")
    print(C)
    print("\nAfter reduction to observer canonical form using duality,")
    print("the new output matrix C_new is:")
    print(np.round(C_new, 4))
    
    print("\nThe new output equations y = C_new * x_new are:")
    for i in range(C_new.shape[0]):
        equation = f"y{i+1} = "
        terms = []
        for j in range(C_new.shape[1]):
            # Format the number to be clean
            num = C_new[i, j]
            if abs(num) > 1e-9: # Only include non-zero terms
                # Format to remove trailing .0
                num_str = f"{num:.4g}".rstrip('0').rstrip('.') if '.' in f"{num:.4g}" else f"{num:.4g}"
                terms.append(f"{num_str}*x_new{j+1}")
        if not terms:
            equation += "0"
        else:
            equation += " + ".join(terms).replace("+ -", "- ")
        print(equation)
