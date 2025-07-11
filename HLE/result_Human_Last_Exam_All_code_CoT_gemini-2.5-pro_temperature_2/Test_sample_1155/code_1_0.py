import numpy as np

def main():
    """
    This program finds the explicit representation of the normal cone T_F^°(x*)
    for a given feasible set F and a point x*.
    """

    # The point of interest
    x_star = np.array([2, 0, -1])

    # The gradients of the active constraints g_1, g_2, g_3 are computed.
    # The constraints g_1, g_2, g_3 were identified as active in the planning phase.
    # g_1(x) = (x_1 - 1)^2 + x_2^2 - 1
    # g_2(x) = (x_1 - 3)^2 + x_2^2 - 1
    # g_3(x) = x_3 + 1
    
    # Gradient of g_1 at x_star
    grad_g1_x_star = np.array([2 * (x_star[0] - 1), 2 * x_star[1], 0])

    # Gradient of g_2 at x_star
    grad_g2_x_star = np.array([2 * (x_star[0] - 3), 2 * x_star[1], 0])

    # Gradient of g_3 at x_star
    grad_g3_x_star = np.array([0, 0, 1])

    active_gradients = [grad_g1_x_star, grad_g2_x_star, grad_g3_x_star]

    # The normal cone T_F^°(x*) is the conic hull of the active gradients.
    # We will construct the string representation of this cone.
    # T_F^°(x*) = { s = mu_1*v1 + mu_2*v2 + ... | mu_i >= 0 }
    
    header = "The normal cone is the set of non-negative linear combinations of the gradients of the active constraints:"
    
    # Format vector components as integers for cleaner output
    vec1_str = f"({int(active_gradients[0][0])}, {int(active_gradients[0][1])}, {int(active_gradients[0][2])})"
    vec2_str = f"({int(active_gradients[1][0])}, {int(active_gradients[1][1])}, {int(active_gradients[1][2])})"
    vec3_str = f"({int(active_gradients[2][0])}, {int(active_gradients[2][1])}, {int(active_gradients[2][2])})"

    equation = (
        f"T_F^°(x*) = {{ s ∈ R^3 | s = "
        f"μ_1 * {vec1_str}^T + "
        f"μ_2 * {vec2_str}^T + "
        f"μ_3 * {vec3_str}^T, "
        f"with μ_1, μ_2, μ_3 ≥ 0 }}"
    )
    
    print(header)
    print(equation)

if __name__ == "__main__":
    main()