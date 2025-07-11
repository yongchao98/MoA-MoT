import numpy as np

# Step 1: Define the original system matrices
A = np.array([[1, 1, 0], 
              [2, 1, 1], 
              [0, 2, 0]])

C = np.array([[0, 1, 0], 
              [1, 1, 0]])

print("Original System Matrices:")
print("A =\n", A)
print("C =\n", C)
print("-" * 30)

# Step 2: Define the dual system using duality principle
# A_d = A^T, B_d = C^T
A_d = A.T
B_d = C.T

print("Dual System Matrices:")
print("A_d = A^T =\n", A_d)
print("B_d = C^T =\n", B_d)
print("-" * 30)

# Step 3: Determine the controllability indices (mu_1, mu_2) of the dual system.
# These will be equal to the observability indices (v_1, v_2) of the original system.
n = A.shape[0]  # Dimension of the system (n=3)
p = C.shape[0]  # Number of outputs (p=2)

# Get the columns of B_d
b1 = B_d[:, 0].reshape(n, 1)
b2 = B_d[:, 1].reshape(n, 1)

# Calculate terms for the controllability matrix W_c = [b1, b2, A_d*b1, A_d*b2, ...]
Ad_b1 = A_d @ b1

# Check for linear independence by forming a matrix and checking its rank/determinant
# We search in the order b1, b2, Ad_b1, ...
# Vector 1: b1
# Vector 2: b2. The matrix [b1, b2] has rank 2.
# Vector 3: Ad_b1. Let's check if {b1, b2, Ad_b1} are linearly independent.
basis_matrix = np.hstack([b1, b2, Ad_b1])
det_basis = np.linalg.det(basis_matrix)

print("Checking for linear independence to find indices...")
print("The first 3 candidate vectors for the basis are b1, b2, A_d*b1:")
print("Basis Matrix =\n", basis_matrix)
print("Determinant =", round(det_basis, 2))
print("Since the determinant is non-zero, the vectors are linearly independent.")
print("\nThe selected vectors are {b1, Ad_b1} (from input 1) and {b2} (from input 2).")

# The number of vectors associated with each input gives the indices
mu1 = 2  # b1, Ad*b1
mu2 = 1  # b2
v1, v2 = mu1, mu2
print(f"The controllability indices of the dual system are (mu_1, mu_2) = ({mu1}, {mu2}).")
print(f"By duality, the observability indices of the original system are (v_1, v_2) = ({v1}, {v2}).")
print("-" * 30)

# Step 4: Construct the new C matrix (C_o) in observer canonical form
# The structure of C_o is determined by the observability indices.
# Row j has a 1 at column v_1 + ... + v_j
C_o = np.zeros((p, n))

# First row has a 1 at column v_1
C_o[0, v1 - 1] = 1

# Second row has a 1 at column v_1 + v_2
C_o[1, v1 + v2 - 1] = 1

print("The new matrix C in observer canonical form (C_o) is:")
print(C_o)
print("-" * 30)

# Step 5: Display the final output equation y_o(k) = C_o * x_o(k)
print("The new output equation y_o(k) = C_o * x_o(k) is given by:")
for i in range(p):
    eq_str = f"y_o_{i+1}(k) = "
    terms = []
    for j in range(n):
        terms.append(f"{int(C_o[i, j])}*x_o_{j+1}(k)")
    eq_str += " + ".join(terms)
    print(eq_str)

<<<[[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]>>>