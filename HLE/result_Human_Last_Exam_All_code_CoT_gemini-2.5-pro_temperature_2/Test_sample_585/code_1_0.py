import numpy as np
from scipy.linalg import expm

def skew_symmetric(v):
    """
    Creates a skew-symmetric matrix from a 3-element vector.
    v: A 3-element numpy array.
    """
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

# Step 1: Define an example attitude deviation vector 'delta' and its covariance 'Sigma'.
# delta is a 3-element vector representing rotation in radians.
delta = np.array([0.1, -0.2, 0.3])

# Sigma is the 3x3 covariance matrix associated with delta.
# It must be symmetric and positive semi-definite.
Sigma = np.array([
    [0.01, 0.001, -0.002],
    [0.001, 0.015, 0.003],
    [-0.002, 0.003, 0.02]
])

# Set print options for better readability
np.set_printoptions(precision=6, suppress=True)

print("This script calculates the post-reset covariance for an attitude deviation in a Kalman filter.\n")
print(f"The initial attitude deviation vector (delta) is:\n{delta}\n")
print(f"The initial covariance matrix (Sigma) is:\n{Sigma}\n")

# Step 2: Compute the skew-symmetric matrix 'delta_hat' from 'delta'.
delta_hat = skew_symmetric(delta)

# Step 3: Compute the rotation matrix 'R_delta' using the matrix exponential.
# This corresponds to the Rodrigues' rotation formula.
R_delta = expm(delta_hat)

print("The transformation is based on the formula: Sigma_prime = R_delta^T * Sigma * R_delta\n")
print(f"First, we compute the rotation matrix R_delta = exp(hat(delta)):\n{R_delta}\n")

# Step 4: Compute the post-reset covariance 'Sigma_prime'.
# The transformation formula is Sigma' = R_delta^T * Sigma * R_delta
Sigma_prime = R_delta.T @ Sigma @ R_delta

# Step 5: Print the final equation with all the numerical values.
print("Now, we compute the post-reset covariance Sigma_prime.\n")
print("The final equation is:\nSigma_prime =\n")
print(f"{R_delta.T}\n\n    * (multiplied by)\n\n{Sigma}\n\n    * (multiplied by)\n\n{R_delta}\n")
print("\nResulting in Sigma_prime =\n")
print(Sigma_prime)

# The final answer is the expression itself
final_expression = "Sigma' = exp(hat(-delta)) * Sigma * exp(hat(delta))"
print(f"\n<<<The exact expression is: {final_expression}>>>")
