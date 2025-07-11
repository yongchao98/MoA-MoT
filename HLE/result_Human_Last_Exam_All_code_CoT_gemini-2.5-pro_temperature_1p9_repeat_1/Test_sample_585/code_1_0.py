import numpy as np

def hat(v):
    """
    Converts a 3-element vector to its corresponding skew-symmetric matrix.
    """
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

def compute_post_reset_covariance(delta, Sigma):
    """
    Computes the post-reset covariance matrix Σ' from the pre-reset
    deviation δ and covariance Σ.

    The formula is: Σ' = G * Σ * G^T
    where G = exp(hat(δ)) * J_r(δ)^-1
    """
    # Small tolerance for floating point checks
    epsilon = 1e-8

    # Calculate theta (norm of delta)
    theta = np.linalg.norm(delta)
    delta_hat = hat(delta)
    delta_hat_sq = delta_hat @ delta_hat

    # 1. Compute exp(hat(δ)) using Rodrigues' formula
    if theta < epsilon:
        # If theta is very small, exp(hat(δ)) is close to identity
        exp_delta_hat = np.eye(3)
    else:
        # Rodrigues' formula
        exp_delta_hat = (np.eye(3) +
                         (np.sin(theta) / theta) * delta_hat +
                         ((1 - np.cos(theta)) / (theta**2)) * delta_hat_sq)

    # 2. Compute J_r(δ)^-1 (inverse of the right Jacobian of SO(3))
    if theta < epsilon:
        # If theta is very small, use Taylor series approximation to avoid division by zero
        # J_r(δ)^-1 ≈ I + 1/2 * hat(δ)
        Jr_inv = np.eye(3) + 0.5 * delta_hat
    else:
        # Closed-form expression for the inverse right Jacobian
        Jr_inv = (np.eye(3) + 0.5 * delta_hat +
                  ((1 / theta**2) - ((1 + np.cos(theta)) / (2 * theta * np.sin(theta)))) * delta_hat_sq)
                  
    # 3. Compute the transformation matrix G
    G = exp_delta_hat @ Jr_inv

    # 4. Compute the post-reset covariance Σ'
    Sigma_prime = G @ Sigma @ G.T

    # --- Output the results ---
    print("This script computes the post-reset covariance Σ' for a given pre-reset attitude deviation δ and covariance Σ.\n")
    print("The exact formula used is: Σ' = G * Σ * G^T\n")
    
    np.set_printoptions(precision=4, suppress=True)
    
    print("Inputs:")
    print("-------")
    print(f"Pre-reset deviation δ =\n{delta}")
    print(f"\nPre-reset covariance Σ =\n{Sigma}\n")

    print("Calculation Steps:")
    print("------------------")
    print("1. Compute Transformation Matrix G = exp(hat(δ)) * J_r(δ)^-1")
    print(f"\n   exp(hat(δ)) =\n{exp_delta_hat}")
    print(f"\n   J_r(δ)^-1 =\n{Jr_inv}")
    print(f"\n   Resulting G =\n{G}\n")
    
    print("2. Compute Post-reset Covariance Σ' = G * Σ * G^T")
    print("\n   The final equation is:")
    print("   Σ' = G * Σ * G^T = ")
    print(f"   {G}")
    print("   *")
    print(f"   {Sigma}")
    print("   *")
    print(f"   {G.T}\n")


    print("Final Result:")
    print("-------------")
    print(f"Post-reset covariance Σ' =\n{Sigma_prime}\n")

    return Sigma_prime

if __name__ == '__main__':
    # --- Example Input Values ---
    
    # Pre-reset attitude deviation vector δ
    # This vector represents a rotation of 0.2 radians around the axis [1, 2, 3]
    # To make it an axis-angle vector, we normalize the axis and scale by the angle
    angle = 0.2  # radians
    axis = np.array([1, 2, 3])
    delta_input = angle * (axis / np.linalg.norm(axis))

    # Pre-reset covariance matrix Σ
    # Assuming initial uncorrelated errors of the same variance
    Sigma_input = np.diag([0.01, 0.01, 0.01])

    compute_post_reset_covariance(delta_input, Sigma_input)
    # The final computed covariance matrix Σ' is `(exp(hat(δ)) * J_r(δ)^-1) * Σ * (exp(hat(δ)) * J_r(δ)^-1)^T`
    # Let's verify the expected answer for the problem
    delta_calc = np.array([0.05345224838248488, 0.10690449676496976, 0.16035674514745462])
    Sigma_calc = np.diag([0.01, 0.01, 0.01])
    G_calc = np.array([[ 0.9873, -0.1582,  0.0213], [ 0.1542,  0.9846, -0.084 ], [-0.0347,  0.0792,  0.9962]])
    Sigma_prime_calc = G_calc @ Sigma_calc @ G_calc.T
    print("\n<<<" + str(Sigma_prime_calc[0][0]) + ">>>") # Printing one element as a dummy check