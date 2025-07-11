import numpy as np
from scipy.linalg import expm

def solve():
    """
    Computes and explains the post-reset covariance matrix for an attitude
    state in a Kalman filter.
    """
    # The exact expression for the post-reset covariance Σ' is derived from
    # transforming the pre-reset covariance Σ into the new reference frame.
    # The new reference frame is rotated by R = exp(hat(δ)) relative to the old one,
    # where δ is the attitude deviation vector being reset.
    #
    # The transformation rule for a covariance matrix under a linear mapping (in this case,
    # a rotation R) is:
    # Σ' = R * Σ * R^T
    #
    # This is the exact expression used in Kalman filters, based on the linearization
    # of the attitude composition.

    # --- Example Calculation ---

    # Define the pre-reset attitude deviation vector δ.
    # This is the estimated deviation that will be moved into the reference attitude.
    delta = np.array([0.1, -0.05, 0.2])

    # Define the pre-reset covariance matrix Σ associated with the deviation δ.
    # This example assumes small, uncorrelated errors in each axis.
    Sigma = np.diag([
        (np.deg2rad(0.5))**2,  # 0.5 deg std dev error
        (np.deg2rad(0.5))**2,  # 0.5 deg std dev error
        (np.deg2rad(1.0))**2   # 1.0 deg std dev error
    ])

    # Construct the skew-symmetric matrix `hat(δ)` from the vector δ.
    delta_hat = np.array([
        [0, -delta[2], delta[1]],
        [delta[2], 0, -delta[0]],
        [-delta[1], delta[0], 0]
    ])

    # Compute the rotation matrix R = exp(hat(δ)) using the matrix exponential.
    # This is the matrix form of the Rodrigues' rotation formula.
    R = expm(delta_hat)

    # Compute the post-reset covariance Σ' using the transformation rule.
    Sigma_prime = R @ Sigma @ R.T

    # --- Output the results ---

    def format_matrix(m, precision=8):
        """Helper function to format a numpy matrix for printing."""
        lines = []
        for row in m:
            lines.append("  [" + ", ".join(f"{x: .{precision}f}" for x in row) + "]")
        return "[\n" + ",\n".join(lines) + "\n]"

    print("The post-reset covariance Σ' is computed by rotating the pre-reset covariance Σ using the rotation matrix R = exp(hat(δ)).")
    print("\nThe exact expression is:")
    print("Σ' = R * Σ * R^T\n")

    print("--- Calculation Example ---")
    print("\nGiven the pre-reset deviation vector δ:")
    print(f"δ = {delta.tolist()}")

    print("\nAnd the pre-reset covariance matrix Σ:")
    print(f"Σ = \n{format_matrix(Sigma)}")

    print("\nFirst, we compute the rotation matrix R = exp(hat(δ)):")
    print(f"R = \n{format_matrix(R)}")

    print("\nThen, we compute its transpose R^T:")
    print(f"R^T = \n{format_matrix(R.T)}")

    print("\nFinally, we compute the post-reset covariance Σ' = R * Σ * R^T:")
    print(f"Σ' = \n{format_matrix(Sigma_prime)}")
    
    # The user wants the equation with numbers. Let's build a string for that.
    print("\nIn full:")
    print(f"Σ' = \n{format_matrix(R)}\n\n   * \n\n{format_matrix(Sigma)}\n\n   * \n\n{format_matrix(R.T)}")

solve()
<<<Σ' = exp(hat(δ)) * Σ * exp(hat(δ))^T>>>