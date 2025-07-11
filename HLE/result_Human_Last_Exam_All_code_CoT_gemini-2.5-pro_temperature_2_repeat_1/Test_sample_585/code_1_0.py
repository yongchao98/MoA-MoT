import numpy as np
from scipy.linalg import expm

def print_equation(sigma_prime, R, sigma):
    """
    Prints the matrices in the equation Σ' = R * Σ * R^T.
    """
    print("Post-reset covariance Σ' is calculated as:\n")
    print("Σ' = R * Σ * R^T\n")
    
    np.set_printoptions(precision=4, suppress=True)

    print("Where R(-δ) =\n", R, "\n")
    print("Σ =\n", sigma, "\n")
    print("R(-δ)^T =\n", R.T, "\n")
    print("Resulting Σ' =\n", sigma_prime, "\n")


def calculate_post_reset_covariance():
    """
    Calculates the post-reset attitude covariance based on the exact expression.
    """
    # 1. Define the pre-reset attitude deviation vector 'delta' that is being
    # absorbed into the reference attitude. This is the filter's current estimate.
    # Example values (in radians)
    delta = np.array([0.1, -0.05, 0.02])

    # 2. Define the corresponding pre-reset covariance matrix 'Sigma'.
    # It must be a 3x3 symmetric positive semi-definite matrix.
    # For this example, we assume some small variances and non-zero covariance.
    Sigma = np.array([[1.5, 0.2, 0.1],
                      [0.2, 1.8, -0.3],
                      [0.1, -0.3, 1.2]]) * 1e-4

    # 3. Compute the skew-symmetric matrix for -delta.
    # The hat operator for a vector v = [vx, vy, vz] gives:
    # [[ 0, -vz,  vy],
    #  [ vz,  0, -vx],
    #  [-vy,  vx,  0]]
    neg_delta_skew = np.array([
        [0,          delta[2], -delta[1]],
        [-delta[2], 0,          delta[0]],
        [delta[1],  -delta[0], 0]
    ])

    # 4. Compute the rotation matrix R(-delta) via the matrix exponential.
    # R(-delta) = expm(-delta_hat)
    R = expm(neg_delta_skew)

    # 5. Compute the post-reset covariance Sigma_prime.
    # Sigma_prime = R * Sigma * R^T
    Sigma_prime = R @ Sigma @ R.T

    # 6. Print the full equation with the numerical matrices
    print_equation(Sigma_prime, R, Sigma)

if __name__ == "__main__":
    calculate_post_reset_covariance()
