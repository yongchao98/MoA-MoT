import numpy as np

def skew(v):
    """
    Converts a 3-element vector to its skew-symmetric matrix form.
    """
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

def compute_post_reset_covariance(delta, Sigma):
    """
    Computes the post-reset covariance matrix using the exact expression.

    Args:
        delta (np.array): The 3-element attitude deviation vector.
        Sigma (np.array): The pre-reset covariance matrix. The first 3 dimensions
                          are assumed to correspond to the attitude deviation.
    """
    print("--- Inputs ---")
    print("Attitude deviation vector (δ):")
    print(delta)
    print("\nPre-reset covariance matrix (Σ):")
    print(Sigma)

    # 1. Calculate the transformation matrix J = exp(skew(δ))^T
    theta = np.linalg.norm(delta)
    if theta < 1e-9:
        # If the rotation is very small, J is the identity matrix
        J = np.identity(3)
    else:
        k = delta / theta
        K = skew(k)
        # Rodrigues' rotation formula for exp(skew(delta))
        R = np.identity(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)
        # The transformation matrix is the transpose
        J = R.T
    
    print("\n--- Calculation ---")
    print("Transformation matrix (J = exp(skew(δ))^T):")
    print(J)

    # 2. Partition the covariance matrix Sigma
    # Assumes attitude is the first 3 states
    Sigma_dd = Sigma[0:3, 0:3]
    Sigma_dr = Sigma[0:3, 3:]
    Sigma_rd = Sigma[3:, 0:3]
    Sigma_rr = Sigma[3:, 3:]

    # 3. Compute the new covariance blocks
    Sigma_dd_new = J @ Sigma_dd @ J.T
    Sigma_dr_new = J @ Sigma_dr
    Sigma_rd_new = Sigma_rd @ J.T
    Sigma_rr_new = Sigma_rr  # Unchanged

    # 4. Assemble the new covariance matrix Sigma'
    Sigma_new = np.block([
        [Sigma_dd_new, Sigma_dr_new],
        [Sigma_rd_new, Sigma_rr_new]
    ])

    print("\n--- Output ---")
    print("Post-reset covariance matrix (Σ'):")
    print(Sigma_new)
    
    # Return the requested single value for the final answer
    return Sigma_new[0, 0]

if __name__ == '__main__':
    # Define an example attitude deviation vector δ
    # These are the small rotation angles in radians to be reset
    delta = np.array([0.1, -0.05, 0.08])

    # Define a sample 6x6 pre-reset covariance matrix Σ
    # This matrix must be symmetric and positive semi-definite.
    # We create a random one for demonstration.
    # The state is [attitude_dev, other_states]
    np.random.seed(42)
    A = np.random.rand(6, 6) * 0.1
    Sigma = A.T @ A + np.diag([0.01, 0.01, 0.01, 0.05, 0.05, 0.05])


    # Compute the post-reset covariance and get the final answer value
    final_answer = compute_post_reset_covariance(delta, Sigma)
    # The final answer format requires a single value. We use Σ'(0,0).
    print(f"\n<<<The value of Σ'(0,0) is: {final_answer}>>>")