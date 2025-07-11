import numpy as np
from scipy.linalg import expm

def solve():
    """
    Calculates the post-reset covariance for an attitude error state in a
    Kalman filter and prints the full equation.
    """
    # This is the estimated state deviation right before the reset step.
    # In the problem description this is simply named δ.
    delta_est = np.array([0.1, -0.1, 0.2])

    # This is the covariance matrix of the error in our estimate of delta.
    # Sigma = Cov(δ_true - δ_est)
    Sigma = np.array([
        [0.01, 0.001, -0.002],
        [0.001, 0.02, 0.005],
        [-0.002, 0.005, 0.03]
    ])

    def skew(v):
        """
        Converts a 3-element vector to its skew-symmetric matrix form.
        """
        return np.array([
            [0, -v[2], v[1]],
            [v[2], 0, -v[0]],
            [-v[1], v[0], 0]
        ])

    # The reset operation updates the reference attitude by R = exp(hat(δ_est)).
    # The covariance must be transformed to be consistent with this new reference.
    # The transformation is Σ' = R^T * Σ * R.

    # 1. Compute the skew-symmetric matrix from the deviation vector.
    Delta_skew = skew(delta_est)

    # 2. Compute the rotation matrix R using the matrix exponential.
    R = expm(Delta_skew)

    # 3. Compute the post-reset covariance Sigma'.
    Sigma_prime = R.T @ Sigma @ R

    # 4. Print the final equation with all numerical values.
    np.set_printoptions(precision=5, suppress=True)

    print("The exact expression for the post-reset covariance Sigma' is:")
    print("Σ' = R^T * Σ * R\n")

    print(f"Given the pre-reset deviation δ_est:\n{delta_est}\n")
    print(f"And the pre-reset covariance Σ:\n{Sigma}\n")

    print(f"First, we compute the rotation matrix R = exp(hat(δ_est)):\nR =\n{R}\n")

    print(f"The transpose of R is R^T:\nR^T =\n{R.T}\n")

    print("Finally, the full calculation for the post-reset covariance Σ' is:\n")
    print("Σ' = \n")
    print(f"{R.T}\n")
    print(" * \n")
    print(f"{Sigma}\n")
    print(" * \n")
    print(f"{R}\n")

    print(f"\nResulting in Σ' =\n{Sigma_prime}")

solve()
<<<Σ' = R^T * Σ * R>>>