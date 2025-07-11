import numpy as np

def compute_post_reset_covariance(delta, Sigma):
    """
    Computes the post-reset attitude error covariance matrix.

    Args:
        delta (np.ndarray): The 3x1 attitude deviation vector.
        Sigma (np.ndarray): The 3x3 pre-reset attitude error covariance matrix.

    Returns:
        np.ndarray: The 3x3 post-reset covariance matrix Sigma_prime.
    """
    # Ensure inputs are numpy arrays
    delta = np.asarray(delta).reshape(3)
    Sigma = np.asarray(Sigma)

    # The norm of the rotation vector
    angle = np.linalg.norm(delta)

    # Handle the case of zero rotation to avoid division by zero
    if angle == 0:
        R = np.identity(3)
    else:
        # The skew-symmetric (cross-product) matrix of delta
        delta_hat = np.array([
            [0, -delta[2], delta[1]],
            [delta[2], 0, -delta[0]],
            [-delta[1], delta[0], 0]
        ])
        
        # Rodrigues' rotation formula: R = I + sin(theta)/theta * K + (1-cos(theta))/theta^2 * K^2
        # where K is the skew-symmetric matrix (delta_hat)
        delta_hat_sq = delta_hat @ delta_hat
        
        # Coefficients for Rodrigues' formula
        a = np.sin(angle) / angle
        b = (1 - np.cos(angle)) / (angle**2)
        
        # The rotation matrix R = exp(delta_hat)
        R = np.identity(3) + a * delta_hat + b * delta_hat_sq

    # The covariance update equation: Sigma' = R * Sigma * R^T
    Sigma_prime = R @ Sigma @ R.T

    # --- Outputting the results ---
    print("The post-reset covariance Sigma' is computed using the formula:")
    print("Sigma' = R * Sigma * R^T\n")
    
    print("Where R is the rotation matrix derived from the deviation delta.")
    print(f"Given delta = {delta.round(4)}\n")
    
    print("The pre-reset covariance matrix (Sigma) is:")
    print(np.round(Sigma, 6))
    print("\n")
    
    print("The computed rotation matrix (R) is:")
    print(np.round(R, 6))
    print("\n")
    
    print("The transpose of the rotation matrix (R^T) is:")
    print(np.round(R.T, 6))
    print("\n")
    
    print("The final post-reset covariance matrix (Sigma') is:")
    print(np.round(Sigma_prime, 6))
    
    return Sigma_prime

if __name__ == '__main__':
    # Example usage:
    # A small attitude deviation vector delta
    delta_k = np.array([0.1, -0.05, 0.02])

    # A sample pre-reset covariance matrix Sigma (assumed diagonal for simplicity)
    Sigma_k = np.diag([0.01**2, 0.01**2, 0.01**2])
    
    # Compute the post-reset covariance
    Sigma_k_prime = compute_post_reset_covariance(delta_k, Sigma_k)
