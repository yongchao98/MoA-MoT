import numpy as np

# Define the system matrices
A = np.array([[1, 1, 0],
              [2, 1, 1],
              [0, 2, 0]])

C = np.array([[0, 1, 0],
              [1, 1, 0]])

# Get the dimension of the state space
n = A.shape[0]

# Get the rows of C, which are the output vectors
c1 = C[0, :]
c2 = C[1, :]

print("Original system matrices:")
print("A =\n", A)
print("C =\n", C)
print("-" * 20)

# Step 1: Determine the observability indices (mu_1, mu_2)
# We check for linear independence of the rows of the observability matrix
# taken in the order: c1, c2, c1*A, c2*A, ...

print("Checking for linear independence to find observability indices:")
# Vector 1: c1
v1 = c1
print("v1 = c1 =", v1)

# Vector 2: c2
# Check if c2 is linearly independent of {v1}
# Since v1 and c2 are not scalar multiples, they are independent.
v2 = c2
print("v2 = c2 =", v2)

# Vector 3: c1*A
# We need to find a third vector to form a basis for R^3
v3 = c1 @ A
print("v3 = c1*A =", v3)

# Check if {v1, v2, v3} are linearly independent by checking the determinant
# of the matrix formed by these vectors.
M = np.vstack([v1, v2, v3])
det_M = np.linalg.det(M)

print("\nMatrix of vectors to test for independence:")
print("M =\n", M)
print("Determinant of M =", det_M)

if abs(det_M) > 1e-9:
    print("\nThe vectors are linearly independent.")
    # We found 3 independent vectors.
    # We used c1 and c1*A from the first output's chain.
    mu1 = 2
    # We used c2 from the second output's chain.
    mu2 = 1
    print(f"The observability indices are mu_1 = {mu1} and mu_2 = {mu2}.")
else:
    print("The vectors are linearly dependent. Stopping.")
    mu1, mu2 = None, None

print("-" * 20)

# Step 2: Construct the new C matrix in observer canonical form
# For a system with observability indices mu_1, mu_2, ...,
# the C matrix in Luenberger observer canonical form is a sparse matrix.
# It has a '1' in the 1st column, the (mu_1 + 1)-th column, etc.
if mu1 is not None:
    # Initialize the new C matrix
    C_new = np.zeros_like(C, dtype=float)

    # For the first output, the new c1_o has a 1 at the 1st position
    C_new[0, 0] = 1

    # For the second output, the new c2_o has a 1 at the (mu1 + 1)-th position
    C_new[1, mu1] = 1

    print("The new matrix C in observer canonical form is:")
    print(C_new)

    print("\nThis is derived from the definition of the Luenberger observer canonical form, where the new output vector y_i selects the first state of the i-th block in the transformed state vector z.")
    print("For y_1, we have y_1 = z_1.")
    print(f"For y_2, the second block starts at index mu_1 + 1 = {mu1+1}, so y_2 = z_{mu1+1}.")
    print("\nFinal equation for the new C matrix:")
    print(f"C_new[0, 0] = {C_new[0, 0]}")
    print(f"C_new[0, 1] = {C_new[0, 1]}")
    print(f"C_new[0, 2] = {C_new[0, 2]}")
    print(f"C_new[1, 0] = {C_new[1, 0]}")
    print(f"C_new[1, 1] = {C_new[1, 1]}")
    print(f"C_new[1, mu1] = C_new[1, 2] = {C_new[1, 2]}")


<<<[[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]>>>