import numpy as np

def hat(v):
    """
    Computes the skew-symmetric matrix (hat operator) for a 3x1 vector.
    """
    return np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]])

def left_jacobian_so3(delta):
    """
    Computes the Left Jacobian of SO(3).
    J_l(d) = I + (1-cos(theta))/theta^2 * [d]x + (theta-sin(theta))/theta^3 * [d]x^2
    """
    delta = np.squeeze(np.asarray(delta))
    theta = np.linalg.norm(delta)
    delta_hat = hat(delta)
    
    if np.isclose(theta, 0.0):
        # Use Taylor series expansion for small theta
        # (1-cos(t))/t^2 -> 1/2
        # (t-sin(t))/t^3 -> 1/6
        A = 0.5
        B = 1.0/6.0
    else:
        theta2 = theta * theta
        theta3 = theta2 * theta
        A = (1 - np.cos(theta)) / theta2
        B = (theta - np.sin(theta)) / theta3

    J_l = np.identity(3) + A * delta_hat + B * (delta_hat @ delta_hat)
    return J_l

# --- Inputs ---
# The attitude deviation vector being reset.
delta = np.array([0.3, -0.2, 0.4])

# The covariance matrix associated with the error in delta.
# This is a symmetric positive semi-definite matrix.
Sigma = np.diag([0.01, 0.01, 0.02])

# --- Calculation ---
# 1. Compute the exact Jacobian for the error transformation
G = left_jacobian_so3(delta)

# 2. Compute the post-reset covariance matrix
Sigma_prime = G @ Sigma @ G.T

# --- Output ---
# Print the final expression with the computed matrices.
# Set print options for better readability
np.set_printoptions(precision=4, suppress=True)

print("The post-reset covariance Σ' is computed as: Σ' = G * Σ * G.T")
print("\nWhere the original covariance Σ is:")
print(Sigma)
print("\nAnd the transformation Jacobian G = J_l(δ) is:")
print(G)
print("\nAnd the resulting post-reset covariance Σ' is:")
print(Sigma_prime)

# Final answer format requirement
print("\n<<<")
print("Σ' = G * Σ * G^T, where G is the left Jacobian J_l(δ). The equation with the given values is:")
print(f"{Sigma_prime} = \n{G} \n* \n{Sigma} \n* \n{G.T}")
print(">>>")