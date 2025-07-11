import math

def display_covariance_reset_formula():
    """
    This function explains and prints the exact formula for the post-reset
    covariance in an error-state Kalman filter for attitude.
    """
    print("The post-reset covariance Σ' is computed by transforming the pre-reset covariance Σ.")
    print("This is done using a Jacobian matrix, G, that represents the exact linearization of the non-linear reset operation on the error state.")
    print("The covariance update follows the standard formula for linear transformations:\n")
    print("  Σ' = G * Σ * G^T\n")
    print("---------------------------------------------------------------------")
    print("DEFINITIONS:\n")
    print("  - Σ: The 3x3 covariance matrix of the attitude deviation vector before the reset.")
    print("  - Σ': The 3x3 covariance matrix after the reset.")
    print("  - δ_hat: The 3-element attitude deviation vector (the filter's estimate) being reset.")
    print("  - θ = ||δ_hat||: The magnitude of the deviation vector, representing the angle of rotation.")
    print("  - I: The 3x3 identity matrix.")
    print("  - [δ_hat]_x: The skew-symmetric matrix form of the vector δ_hat.")
    print("    If δ_hat = [dx, dy, dz]^T, then:")
    print("          [  0  -dz   dy  ]")
    print("    [δ_hat]_x = [  dz   0  -dx  ]")
    print("          [ -dy   dx   0   ]\n")
    print("---------------------------------------------------------------------")
    print("EXACT TRANSFORMATION MATRIX (G):\n")
    print("The exact, closed-form expression for the Jacobian matrix G is:\n")
    print("  G = I + (1/2) * [δ_hat]_x + C(θ) * ([δ_hat]_x)^2\n")
    print("where C(θ) is a scalar coefficient that depends on the magnitude θ:\n")
    print("               1          (1 + cos(θ))")
    print("  C(θ) = ---------- - ---------------------")
    print("            θ^2          2 * θ * sin(θ)\n")
    print("Note: For implementation, a Taylor series expansion or L'Hopital's rule must be used for C(θ) when θ is very close to zero to avoid numerical instability.\n")
    print("---------------------------------------------------------------------")
    print("FINAL EQUATION FOR POST-RESET COVARIANCE (Σ'):\n")
    print("Putting it all together, the exact expression is:\n")
    print("  Σ' = (I + (1/2)*[δ_hat]_x + C(θ)*([δ_hat]_x)^2) * Σ * (I + (1/2)*[δ_hat]_x + C(θ)*([δ_hat]_x)^2)^T")

if __name__ == '__main__':
    display_covariance_reset_formula()
