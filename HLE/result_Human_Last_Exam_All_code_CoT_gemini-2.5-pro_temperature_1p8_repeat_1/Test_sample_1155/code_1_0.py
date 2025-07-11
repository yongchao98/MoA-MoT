import numpy as np

def solve_normal_cone():
    """
    This script calculates the explicit representation of the normal cone for a given
    feasible set at a specific point.
    """
    # Step 1: Define the point x* and constraint functions
    x_star = np.array([2, 0, -1])

    g_functions = [
        lambda x: (x[0] - 1)**2 + x[1]**2 - 1,
        lambda x: (x[0] - 3)**2 + x[1]**2 - 1,
        lambda x: x[2] + 1,
        lambda x: -x[2] - 2
    ]
    
    # Gradient functions (Jacobian of g)
    grad_g_functions = [
        lambda x: np.array([2 * (x[0] - 1), 2 * x[1], 0]),
        lambda x: np.array([2 * (x[0] - 3), 2 * x[1], 0]),
        lambda x: np.array([0, 0, 1]),
        lambda x: np.array([0, 0, -1])
    ]

    print(f"The given point is x* = {x_star}\n")

    # Step 2: Identify active constraints
    print("--- Identifying Active Constraints ---")
    g_values = np.array([g(x_star) for g in g_functions])
    
    # Use a small tolerance for floating-point comparison
    tolerance = 1e-9
    active_indices = np.where(np.abs(g_values) < tolerance)[0]

    for i, val in enumerate(g_values):
        is_active = " (Active)" if i in active_indices else ""
        print(f"g_{i+1}(x*) = {val:.4f}{is_active}")
    
    print(f"\nThe active constraint indices are I(x*) = {{ {', '.join(map(str, active_indices + 1))} }}\n")

    # Step 3: Compute gradients of active constraints
    print("--- Computing Gradients of Active Constraints ---")
    active_gradients = [grad_g_functions[i](x_star) for i in active_indices]
    
    for i, grad in enumerate(active_gradients):
        idx = active_indices[i]
        print(f"∇g_{idx+1}(x*) = {grad}")
    
    print("\nSince the feasible set F is convex, the normal cone T_F^°(x^*) is the conic hull of these gradients.\n")

    # Step 4: Describe the normal cone and simplify
    print("--- Explicit Representation of the Normal Cone T_F^°(x^*) ---")
    
    s_eq_parts = []
    for i, grad in enumerate(active_gradients):
        mu_index = i + 1
        s_eq_parts.append(f"μ_{mu_index} * {grad.tolist()}")

    print("A vector s belongs to the normal cone if it can be written as a non-negative linear combination of the active gradients:")
    print(f"s = {' + '.join(s_eq_parts)}, where μ_i ≥ 0 for all i.")
    
    # Constructing the equations for s = (s1, s2, s3)
    print("\nLet s = (s1, s2, s3). This vector equation yields a system of linear equations for its components:")
    
    # Equation for s1
    s1_eq = "s1 = " + " + ".join([f"{grad[0]:.0f}*μ_{active_indices[i]+1}" for i, grad in enumerate(active_gradients) if grad[0] != 0]).replace('+ -', '- ')
    # Equation for s2
    s2_eq = "s2 = " + " + ".join([f"{grad[1]:.0f}*μ_{active_indices[i]+1}" for i, grad in enumerate(active_gradients) if grad[1] != 0])
    # Equation for s3
    s3_eq = "s3 = " + " + ".join([f"{grad[2]:.0f}*μ_{active_indices[i]+1}" for i, grad in enumerate(active_gradients) if grad[2] != 0])
    
    # Manually simplifying based on computed gradients for better output
    grad1_val = active_gradients[0][0]
    grad2_val = active_gradients[1][0]
    grad3_val = active_gradients[2][2]
    print(f"s1 = {grad1_val:.0f}*μ_1 + {grad2_val:.0f}*μ_2")
    print("s2 = 0")
    print(f"s3 = {grad3_val:.0f}*μ_3")

    print("\nSimplifying these conditions:")
    print(f"1. From s3 = μ_3 and μ_3 ≥ 0, we must have s3 ≥ 0.")
    print(f"2. The second component, s2, is always 0.")
    print(f"3. From s1 = {grad1_val:.0f}*μ_1 - {-grad2_val:.0f}*μ_2, and since μ_1, μ_2 can be any non-negative reals, s1 can be any real number.")

    print("\nTherefore, the explicit representation of the normal cone is:")
    print("T_F^°(x^*) = { s = (s1, s2, s3) ∈ R^3 | s2 = 0, s3 ≥ 0 }")

if __name__ == '__main__':
    solve_normal_cone()
    # The final mathematical expression cannot be assigned to a variable, so we print it
    # as the conclusion of our derivation.
    final_answer = "{ s = (s1, s2, s3) in R^3 | s2 = 0, s3 >= 0 }"
    # print(f"\n<<<{final_answer}>>>") # This would be the final answer format