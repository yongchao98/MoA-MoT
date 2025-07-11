import numpy as np

def skew(v):
    """
    Converts a 3-element vector to its skew-symmetric matrix.
    """
    if not hasattr(v, "__len__") or len(v) != 3:
        raise ValueError("Input must be a 3-element vector.")
    return np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]], dtype=float)

def rodrigues_exp(v):
    """
    Computes the matrix exponential of a skew-symmetric matrix using Rodrigues' formula.
    exp(hat(v)) = I + (sin(theta)/theta)*hat(v) + ((1-cos(theta))/theta^2)*hat(v)^2
    """
    theta = np.linalg.norm(v)
    if theta < 1e-12: # Handle the case of a very small rotation
        # Use Taylor series expansion for small theta
        # sin(theta)/theta approx 1 - theta^2/6
        # (1-cos(theta))/theta^2 approx 1/2 - theta^2/24
        A = 1.0 - theta**2 / 6.0
        B = 0.5 - theta**2 / 24.0
    else:
        A = np.sin(theta) / theta
        B = (1.0 - np.cos(theta)) / (theta**2)

    v_skew = skew(v)
    I = np.identity(3)
    R = I + A * v_skew + B * np.dot(v_skew, v_skew)
    return R

def compute_post_reset_covariance(delta, Sigma):
    """
    Computes the post-reset covariance matrix Σ'.
    
    Args:
        delta (np.array): The 3D attitude deviation vector used for the reset.
        Sigma (np.array): The 3x3 covariance matrix of delta before the reset.
    """
    print("--- Input Data ---")
    print(f"Attitude deviation vector δ:\n{delta}\n")
    print(f"Pre-reset covariance Σ:\n{Sigma}\n")
    
    # The transformation matrix is G = exp(-hat(δ))
    # Note: exp(-hat(δ)) is the same as exp(hat(-δ))
    G = rodrigues_exp(-delta)
    G_T = G.T
    
    # The post-reset covariance is Σ' = G * Σ * G^T
    Sigma_prime = G @ Sigma @ G_T
    
    print("--- Calculation ---")
    print("The post-reset covariance Σ' is calculated as:")
    print("Σ' = G * Σ * G^T\n")
    print("where G = exp(-hat(δ)) is the transformation matrix.\n")
    
    np.set_printoptions(precision=4, suppress=True)
    
    print(f"Computed Transformation Matrix G = exp(-hat(δ)):\n{G}\n")
    
    print("--- Final Equation with Numbers ---")
    print("Σ' = \n")
    print(f"{G}\n")
    print("     * \n")
    print(f"{Sigma}\n")
    print("     * \n")
    print(f"{G_T}\n")
    
    print("\n--- Final Result ---")
    print(f"Post-reset covariance Σ':\n{Sigma_prime}\n")

if __name__ == '__main__':
    # Example usage:
    # Let's assume the filter estimated this deviation vector
    delta_reset = np.array([0.1, -0.05, 0.02])

    # And this is the associated 3x3 covariance matrix from the full state covariance
    # (assuming some small, non-diagonal covariance)
    Sigma_pre_reset = np.array([[0.0025, 0.0001, -0.0002],
                                [0.0001, 0.0030,  0.0003],
                                [-0.0002, 0.0003, 0.0020]])

    compute_post_reset_covariance(delta_reset, Sigma_pre_reset)