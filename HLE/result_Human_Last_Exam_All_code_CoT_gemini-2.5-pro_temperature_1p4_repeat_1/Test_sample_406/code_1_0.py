import sys

def check_fgh_tripled_fixed_point_conditions(alpha, beta, gamma):
    """
    Checks the contractive condition for an FGH-tripled fixed point theorem.

    This theorem provides conditions for the existence and uniqueness of a point (x, y, z)
    in complete metric spaces (X, Y, Z) such that:
    F(x, y, z) = x
    G(y, x, y) = y
    H(z, y, x) = z

    The core condition is that the functions F, G, H must be contractive
    with constants alpha, beta, gamma, whose grouped sum is less than 1.

    Args:
        alpha (list or tuple): A list of 3 floats [α₁, α₂, α₃] >= 0.
        beta (list or tuple): A list of 3 floats [β₁, β₂, β₃] >= 0.
        gamma (list or tuple): A list of 3 floats [γ₁, γ₂, γ₃] >= 0.
    """
    if not all(c >= 0 for c in alpha + beta + gamma):
        print("Error: All contractive constants must be non-negative.")
        return

    # Unpack constants for clarity
    alpha_1, alpha_2, alpha_3 = alpha
    beta_1, beta_2, beta_3 = beta
    gamma_1, gamma_2, gamma_3 = gamma

    # The theorem requires a contraction on the product space X x Y x Z.
    # We define an operator T(x,y,z) = (F(x,y,z), G(y,x,y), H(z,y,x)).
    # The condition for T to be a contraction is that max(k_x, k_y, k_z) < 1.
    
    # Calculate k_x, the coefficient for the distance in X
    print("Calculating k_x...")
    print(f"k_x = α₁ + β₂ + γ₃")
    k_x = alpha_1 + beta_2 + gamma_3
    print(f"k_x = {alpha_1} + {beta_2} + {gamma_3} = {k_x:.4f}\n")

    # Calculate k_y, the coefficient for the distance in Y
    print("Calculating k_y...")
    print(f"k_y = α₂ + β₁ + β₃ + γ₂")
    k_y = alpha_2 + beta_1 + beta_3 + gamma_2
    print(f"k_y = {alpha_2} + {beta_1} + {beta_3} + {gamma_2} = {k_y:.4f}\n")

    # Calculate k_z, the coefficient for the distance in Z
    print("Calculating k_z...")
    print(f"k_z = α₃ + γ₁")
    k_z = alpha_3 + gamma_1
    print(f"k_z = {alpha_3} + {gamma_1} = {k_z:.4f}\n")
    
    # The overall contraction constant k is the maximum of these
    k = max(k_x, k_y, k_z)
    
    print("Final Check:")
    print(f"k = max(k_x, k_y, k_z) = max({k_x:.4f}, {k_y:.4f}, {k_z:.4f}) = {k:.4f}")

    if k < 1:
        print(f"Result: Condition met (k = {k:.4f} < 1).")
        print("A unique FGH-tripled fixed point exists, provided the spaces are complete.")
    else:
        print(f"Result: Condition NOT met (k = {k:.4f} >= 1).")
        print("This theorem does not guarantee the existence of a unique FGH-tripled fixed point.")


if __name__ == '__main__':
    # Example 1: A set of constants that satisfy the condition
    print("--- Example 1: Condition Met ---")
    alpha_coeffs_1 = [0.2, 0.1, 0.1]
    beta_coeffs_1 = [0.15, 0.2, 0.1]
    gamma_coeffs_1 = [0.2, 0.1, 0.25]
    check_fgh_tripled_fixed_point_conditions(alpha_coeffs_1, beta_coeffs_1, gamma_coeffs_1)
    
    print("\n" + "="*40 + "\n")
    
    # Example 2: A set of constants that do NOT satisfy the condition
    print("--- Example 2: Condition Not Met ---")
    alpha_coeffs_2 = [0.4, 0.3, 0.2]
    beta_coeffs_2 = [0.2, 0.4, 0.2]
    gamma_coeffs_2 = [0.3, 0.5, 0.4]
    check_fgh_tripled_fixed_point_conditions(alpha_coeffs_2, beta_coeffs_2, gamma_coeffs_2)