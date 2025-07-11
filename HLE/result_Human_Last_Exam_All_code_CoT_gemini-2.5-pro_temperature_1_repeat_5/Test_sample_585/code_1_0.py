import numpy as np
from scipy.linalg import expm

def solve():
    """
    Calculates the post-reset covariance for an attitude deviation in a Kalman filter.
    """
    print("This script calculates the post-reset covariance Σ' based on the formula:")
    print("Σ' = R * Σ * R^T\n")
    print("where R = exp(hat(δ)) is the rotation matrix for the reset, and Σ is the pre-reset covariance.\n")

    # Define an example pre-reset attitude deviation vector δ.
    # This vector represents the mean of the error distribution to be reset.
    delta = np.array([0.1, -0.05, 0.02])

    # Define an example pre-reset covariance matrix Σ.
    # This represents the uncertainty of the attitude deviation.
    # (e.g., corresponding to 1-sigma uncertainties of 1.0, 1.2, and 1.5 degrees)
    sigmas_rad = np.deg2rad([1.0, 1.2, 1.5])
    Sigma = np.diag(sigmas_rad**2)

    # --- Calculation ---

    # 1. Form the skew-symmetric matrix `hat(δ)` from the vector `δ`
    delta_skew = np.array([
        [0,        -delta[2],  delta[1]],
        [delta[2],  0,         -delta[0]],
        [-delta[1], delta[0],  0]
    ])

    # 2. Compute the rotation matrix R = exp(hat(δ)) using the matrix exponential.
    #    scipy.linalg.expm is a robust method for this.
    R = expm(delta_skew)

    # 3. Compute the post-reset covariance Σ' using the derived formula.
    #    R.T is the transpose of R, which is equal to exp(-hat(δ)).
    Sigma_prime = R @ Sigma @ R.T

    # --- Output the results ---
    print("--- Example Calculation ---\n")
    print(f"Given the attitude deviation vector δ to be reset:\n{delta}\n")
    print(f"And its pre-reset covariance Σ:\n{Sigma}\n")
    print(f"The corresponding rotation matrix R = exp(hat(δ)) is:\n{R}\n")
    print("The final equation is Σ' = R * Σ * R.T")
    print(f"The resulting post-reset covariance Σ' is:\n{Sigma_prime}\n")

solve()
<<<Σ' = exp(hat(δ)) * Σ * exp(-hat(δ)))>>>