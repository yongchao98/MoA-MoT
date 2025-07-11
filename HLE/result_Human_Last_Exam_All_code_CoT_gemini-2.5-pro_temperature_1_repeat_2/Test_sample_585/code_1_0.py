import numpy as np

def hat(v):
    """
    Computes the skew-symmetric matrix (hat map) for a 3D vector.
    """
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

def J_l(phi):
    """
    Computes the exact Left Jacobian of SO(3).
    
    Args:
        phi: A 3D rotation vector.
        
    Returns:
        The 3x3 Left Jacobian matrix.
    """
    theta = np.linalg.norm(phi)
    
    # Handle the case of a very small rotation to avoid division by zero
    if theta < 1e-9:
        return np.identity(3)
        
    phi_x = hat(phi)
    phi_x2 = phi_x @ phi_x
    
    # Coefficients from the closed-form expression
    A = (1 - np.cos(theta)) / (theta**2)
    B = (theta - np.sin(theta)) / (theta**3)
    
    return np.identity(3) + A * phi_x + B * phi_x2

def main():
    """
    Calculates and prints the post-reset covariance matrix.
    """
    # --- Inputs ---
    # δ: The mean of the attitude deviation vector before the reset.
    # Let's assume a rotation of roughly 0.37 radians (21 degrees).
    delta = np.array([0.1, -0.2, 0.3])

    # Σ: The covariance matrix associated with δ.
    # Let's assume a standard deviation of 0.05 radians (about 3 degrees) on each axis.
    stdev = 0.05
    Sigma = np.diag([stdev**2, stdev**2, stdev**2])

    # --- Calculation ---
    # 1. Compute the Left Jacobian for the given delta
    Jl = J_l(delta)

    # 2. Compute the post-reset covariance using the exact formula
    Sigma_prime = Jl @ Sigma @ Jl.T

    # --- Output ---
    np.set_printoptions(precision=6, suppress=True)
    
    print("The post-reset covariance Σ' is computed using the formula:")
    print("Σ' = J_l(δ) * Σ * J_l(δ)ᵀ\n")
    print("--- Each Number in the Final Equation ---")
    
    print("δ (pre-reset attitude deviation mean):")
    print(delta)
    
    print("\nJ_l(δ) (Left Jacobian of SO(3) for δ):")
    print(Jl)

    print("\nΣ (pre-reset covariance):")
    print(Sigma)
    
    print("\nJ_l(δ)ᵀ (Transpose of the Jacobian):")
    print(Jl.T)

    print("\nResulting Σ' (post-reset covariance):")
    print(Sigma_prime)
    
    # For a direct answer format, we can print the final matrix.
    # Convert the matrix to a string to fit the requested format.
    final_answer_str = np.array2string(Sigma_prime, separator=',', formatter={'float_kind':lambda x: "%.6f" % x})
    # This is just for demonstration, the primary output is above.
    # print(f"\n<<<Final Answer Matrix:\n{final_answer_str}>>>")


if __name__ == "__main__":
    main()