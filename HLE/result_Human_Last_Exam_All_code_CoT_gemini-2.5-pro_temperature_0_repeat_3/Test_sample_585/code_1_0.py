import numpy as np
from scipy.linalg import expm

def calculate_post_reset_covariance(delta, Sigma):
    """
    Calculates the post-reset covariance matrix for an attitude deviation.

    The reset operation transforms the covariance according to the exact formula:
    Sigma_prime = exp(hat(-delta/2)) * Sigma * exp(hat(delta/2))

    Args:
        delta (np.ndarray): The 3x1 attitude deviation vector being reset.
        Sigma (np.ndarray): The 3x3 pre-reset covariance matrix of the deviation.

    Returns:
        np.ndarray: The 3x3 post-reset covariance matrix.
    """
    # Ensure inputs are numpy arrays
    delta = np.asarray(delta)
    Sigma = np.asarray(Sigma)

    if delta.shape != (3,) or Sigma.shape != (3, 3):
        raise ValueError("Invalid input shapes. delta must be (3,) and Sigma must be (3,3).")

    def skew(v):
        """
        Converts a 3-element vector to its skew-symmetric matrix form (hat map).
        """
        return np.array([
            [0, -v[2], v[1]],
            [v[2], 0, -v[0]],
            [-v[1], v[0], 0]
        ])

    # The formula requires rotating by -delta/2 and delta/2
    v1 = -delta / 2.0
    v2 = delta / 2.0

    # Get the skew-symmetric matrices
    v1_skew = skew(v1)
    v2_skew = skew(v2)

    # The transformation matrices are the matrix exponentials of the skew-symmetric matrices.
    # scipy.linalg.expm provides a robust implementation of the matrix exponential.
    G = expm(v1_skew)
    G_T = expm(v2_skew)

    # Apply the covariance transformation law: Sigma' = G * Sigma * G^T
    Sigma_prime = G @ Sigma @ G_T

    # --- Output the results ---
    print("The exact expression for the post-reset covariance Σ' is:")
    print("Σ' = exp(hat(-δ/2)) * Σ * exp(hat(δ/2))\n")

    print("Given the pre-reset deviation vector:")
    print(f"δ = \n{delta}\n")

    print("And the pre-reset covariance matrix:")
    print(f"Σ = \n{Sigma}\n")

    print("The first transformation matrix is G = exp(hat(-δ/2)):")
    print(f"G = \n{G}\n")

    print("The second transformation matrix is G_transpose = exp(hat(δ/2)):")
    print(f"G_transpose = \n{G_T}\n")

    print("The final post-reset covariance matrix is Σ' = G * Σ * G_transpose:")
    print(f"Σ' = \n{Sigma_prime}\n")

    return Sigma_prime

if __name__ == '__main__':
    # Example Usage:
    # Define a sample pre-reset attitude deviation vector (in radians)
    # This represents a small rotation that the filter has estimated.
    delta_k = np.array([0.08, -0.05, 0.12])

    # Define a sample pre-reset covariance matrix for this deviation.
    # This typically has small values on the diagonal, representing low uncertainty.
    Sigma_k = np.diag([np.deg2rad(0.1)**2, np.deg2rad(0.1)**2, np.deg2rad(0.15)**2])

    # Calculate and print the post-reset covariance
    calculate_post_reset_covariance(delta_k, Sigma_k)