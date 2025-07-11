import numpy as np

def skew_symmetric(v):
    """
    Computes the skew-symmetric matrix for a 3-element vector v.
    hat(v) = [[ 0, -v2,  v1],
              [ v2,   0, -v0],
              [-v1,  v0,   0]]
    """
    return np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]])

def compute_post_reset_covariance(delta, Sigma):
    """
    Computes the post-reset covariance matrix for an attitude deviation.

    Args:
        delta (np.ndarray): The 3-element attitude deviation vector (δ).
        Sigma (np.ndarray): The 3x3 pre-reset covariance matrix (Σ).

    Returns:
        np.ndarray: The 3x3 post-reset covariance matrix (Σ').
    """
    # Ensure inputs are numpy arrays
    delta = np.asarray(delta)
    Sigma = np.asarray(Sigma)

    if delta.shape != (3,) or Sigma.shape != (3, 3):
        raise ValueError("Inputs must be a 3-element vector and a 3x3 matrix.")

    # Calculate the norm of delta
    theta = np.linalg.norm(delta)

    # Compute the skew-symmetric matrix of delta
    delta_skew = skew_symmetric(delta)

    # Compute the transformation matrix G = exp(-hat(delta)) using Rodrigues' formula
    if theta < 1e-9:  # Handle the case of a very small rotation
        G = np.identity(3)
    else:
        I = np.identity(3)
        # Note: Rodrigues' formula for exp(-hat(delta)) has a minus sign on the sin term
        # G = I - (sin(theta)/theta) * hat(delta) + ((1-cos(theta))/theta^2) * hat(delta)^2
        c1 = np.sin(theta) / theta
        c2 = (1 - np.cos(theta)) / (theta**2)
        G = I - c1 * delta_skew + c2 * (delta_skew @ delta_skew)

    # Compute the post-reset covariance: Sigma' = G * Sigma * G^T
    Sigma_prime = G @ Sigma @ G.T
    
    return G, Sigma_prime

# --- Example Usage ---
# Define an example attitude deviation vector delta
delta_k = np.array([0.1, -0.2, 0.15])

# Define an example pre-reset covariance matrix Sigma
# Using a non-identity matrix to show a more general case
Sigma_k = np.array([[1.0, 0.1, 0.2],
                    [0.1, 2.0, 0.3],
                    [0.2, 0.3, 3.0]])

# Compute the transformation matrix G and the post-reset covariance Sigma'
G_matrix, Sigma_k_prime = compute_post_reset_covariance(delta_k, Sigma_k)

# --- Print the results in equation form ---
np.set_printoptions(precision=4, suppress=True)

print("The post-reset covariance Σ' is computed using the exact formula:")
print("Σ' = G * Σ * G^T\n")

print("Where:")
print("δ (attitude deviation vector) =\n", delta_k, "\n")
print("Σ (pre-reset covariance) =\n", Sigma_k, "\n")
print("G = exp(-hat(δ)) (transformation matrix) =\n", G_matrix, "\n")
print("The resulting post-reset covariance Σ' is:\n")
print(Sigma_k_prime)
