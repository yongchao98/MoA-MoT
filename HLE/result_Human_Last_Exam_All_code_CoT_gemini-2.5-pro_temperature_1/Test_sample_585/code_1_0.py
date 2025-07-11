import numpy as np

def skew(v):
    """
    Converts a 3-element vector to its skew-symmetric matrix form.
    [v]ₓ
    """
    if not isinstance(v, np.ndarray) or v.shape != (3,):
        raise ValueError("Input must be a 3-element NumPy array.")
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

def left_jacobian_so3(delta):
    """
    Computes the left Jacobian of SO(3) for a given rotation vector delta.
    This is the exact expression, J_l(delta).
    """
    if not isinstance(delta, np.ndarray) or delta.shape != (3,):
        raise ValueError("Input must be a 3-element NumPy array.")
        
    theta = np.linalg.norm(delta)
    delta_skew = skew(delta)
    I = np.identity(3)

    # Use Taylor series expansion for small theta to avoid division by zero
    # and maintain numerical stability.
    if theta < 1e-9:
        # For theta -> 0, J_l -> I + 1/2*[delta]x
        # For even smaller theta, J_l -> I
        return I + 0.5 * delta_skew
    else:
        theta2 = theta * theta
        theta3 = theta2 * theta
        
        A = (1 - np.cos(theta)) / theta2
        B = (theta - np.sin(theta)) / theta3
        
        # The exact formula for the left Jacobian
        J_l = I + A * delta_skew + B * (delta_skew @ delta_skew)
        return J_l

# --- Main execution ---
# Define a sample attitude deviation vector `delta` (rad) to be reset.
# This would be the filter's state estimate before the reset.
delta = np.array([0.1, -0.2, 0.15])

# Define a sample pre-reset covariance matrix `Sigma` for the deviation `delta`.
# This would be the corresponding block from the filter's covariance matrix.
# For this example, we assume small, uncorrelated errors in each axis.
Sigma = np.diag([0.01**2, 0.01**2, 0.01**2])

# 1. Calculate the exact left Jacobian J_l(delta)
J_l = left_jacobian_so3(delta)

# 2. Compute the post-reset covariance Sigma_prime using the exact formula
# Sigma' = J_l(delta) * Sigma * J_l(delta)^T
Sigma_prime = J_l @ Sigma @ J_l.T

# --- Output the results ---
print("The post-reset covariance Σ' is computed using the formula:")
print("Σ' = J_l(δ) * Σ * J_l(δ)ᵀ\n")

print("Given the pre-reset deviation vector δ:")
print(delta)
print("\nAnd the pre-reset covariance matrix Σ:")
print(Sigma)

print("\nFirst, we compute the left Jacobian of SO(3), J_l(δ):")
print(J_l)

print("\nFinally, we compute the post-reset covariance matrix Σ':")
print(Sigma_prime)

print("\n<<<Σ' = J_l(δ) * Σ * J_l(δ)ᵀ>>>")