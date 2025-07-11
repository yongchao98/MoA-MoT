import numpy as np

def skew(v):
    """
    Converts a 3-element vector to its 3x3 skew-symmetric matrix.
    
    Args:
        v (np.ndarray): A 3-element vector.
        
    Returns:
        np.ndarray: The 3x3 skew-symmetric matrix.
    """
    # Ensure v is a 1D array
    v = np.asarray(v).flatten()
    return np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]])

def compute_post_reset_covariance(delta, Sigma):
    """
    Computes the post-reset covariance matrix for an attitude error state
    using the exact Jacobian transformation.

    Args:
        delta (np.ndarray): The 3x1 attitude deviation vector being reset.
        Sigma (np.ndarray): The 3x3 pre-reset covariance matrix.

    Returns:
        np.ndarray: The 3x3 post-reset covariance matrix.
    """
    # Ensure inputs are numpy arrays
    delta = np.asarray(delta).reshape(3, 1)
    Sigma = np.asarray(Sigma)
    
    # Calculate the norm of the rotation vector
    theta = np.linalg.norm(delta)
    
    # Get the skew-symmetric matrix and its square
    delta_skew = skew(delta)
    delta_skew_sq = delta_skew @ delta_skew

    # Calculate the right Jacobian Jr using Rodrigues' formula.
    # A numerically stable implementation is used for small theta.
    if theta < 1e-6:
        # Use Taylor series expansion for small theta to avoid division by zero.
        # Jr ≈ I - 1/2*δ̂ + 1/6*δ̂²
        A = 1/2 - theta**2 / 24
        B = 1/6 - theta**2 / 120
        Jr = np.identity(3) - A * delta_skew + B * delta_skew_sq
    else:
        theta2 = theta * theta
        theta3 = theta2 * theta
        c_theta = np.cos(theta)
        s_theta = np.sin(theta)
        
        # Coefficients from the Rodrigues' Jacobian formula
        A = (1 - c_theta) / theta2
        B = (theta - s_theta) / theta3
        
        Jr = np.identity(3) - A * delta_skew + B * delta_skew_sq

    # Propagate the covariance: Σ' = J * Σ * Jᵀ
    Sigma_new = Jr @ Sigma @ Jr.T
    
    print("This script computes the post-reset covariance matrix Σ' using the formula:")
    print("Σ' = J * Σ * Jᵀ")
    print("\nWhere J is the right Jacobian of SO(3), J_r(δ).")
    
    print("\n--- Input Values ---")
    print("Reset vector δ:")
    print(delta.flatten())
    print("\nPre-reset covariance Σ:")
    print(Sigma)
    
    print("\n--- Calculated Matrices (The Equation) ---")
    print("J = J_r(δ):")
    print(Jr)
    print("\nΣ:")
    print(Sigma)
    print("\nJᵀ:")
    print(Jr.T)
    
    print("\n--- Final Result ---")
    print("Post-reset covariance Σ':")
    print(Sigma_new)
    
    return Sigma_new

# --- Example Usage ---
# Define an example attitude deviation vector δ_hat to be reset
delta_hat = np.array([0.1, -0.2, 0.15])

# Define an example pre-reset covariance matrix Σ
# (e.g., uncertainties of ~0.5 deg, 0.5 deg, 1.0 deg)
stdev = np.deg2rad([0.5, 0.5, 1.0])
Sigma_pre_reset = np.diag(stdev**2)

# Compute and print the post-reset covariance and all components
Sigma_post_reset = compute_post_reset_covariance(delta_hat, Sigma_pre_reset)
