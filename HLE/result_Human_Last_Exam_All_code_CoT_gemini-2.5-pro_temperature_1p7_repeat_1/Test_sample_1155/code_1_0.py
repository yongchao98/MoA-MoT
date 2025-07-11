import numpy as np

def solve_normal_cone():
    """
    Finds the explicit representation of the normal cone T_F^°(x^*) for a given
    feasible set F and point x^*.
    """
    # Define the point and constraints
    x_star = np.array([2., 0., -1.])

    g_str = [
        "(x_1 - 1)^2 + x_2^2 - 1",
        "(x_1 - 3)^2 + x_2^2 - 1",
        "x_3 + 1",
        "-x_3 - 2"
    ]
    
    g_funcs = [
        lambda x: (x[0] - 1)**2 + x[1]**2 - 1,
        lambda x: (x[0] - 3)**2 + x[1]**2 - 1,
        lambda x: x[2] + 1,
        lambda x: -x[2] - 2
    ]
    
    nabla_g_funcs = [
        lambda x: np.array([2 * (x[0] - 1), 2 * x[1], 0]),
        lambda x: np.array([2 * (x[0] - 3), 2 * x[1], 0]),
        lambda x: np.array([0., 0., 1.]),
        lambda x: np.array([0., 0., -1.])
    ]

    print("--- Problem Setup ---")
    print(f"The point is x* = {x_star.T}")
    print("\nThe inequality constraints g_i(x) <= 0 are:")
    for i, s in enumerate(g_str):
        print(f"g_{i+1}(x) = {s}")
    print("-" * 25 + "\n")

    # Step 1: Identify active constraints at x*
    print("--- Step 1: Identify Active Constraints ---")
    active_indices = []
    print(f"Evaluating constraints at x* = {x_star.T}:")
    for i, g_i in enumerate(g_funcs):
        val = g_i(x_star)
        # A constraint is active if g_i(x*) = 0
        if np.isclose(val, 0):
            status = "-> Active"
            active_indices.append(i)
        else:
            status = "-> Inactive"
        print(f"g_{i+1}(x*) = {val:.4f}  {status}")
    
    print(f"\nThe set of indices of active constraints is I(x*) = {{ {', '.join(str(i+1) for i in active_indices)} }}")
    print("-" * 25 + "\n")

    # Step 2: Compute gradients of active constraints
    print("--- Step 2: Compute Gradients of Active Constraints ---")
    active_gradients = []
    print("Evaluating the gradients of active constraints at x*:")
    for i in active_indices:
        grad = nabla_g_funcs[i](x_star)
        active_gradients.append(grad)
        print(f"nabla g_{i+1}(x*) = {grad.T}")
    print("-" * 25 + "\n")

    # Step 3 & 4: Define and Characterize the Normal Cone
    print("--- Step 3: Define and Characterize the Normal Cone ---")
    print("Since all constraint functions g_i(x) are convex, the feasible set F is convex.")
    print("Therefore, the normal cone is the conic hull of the gradients of the active constraints:")
    print("T_F^°(x^*) = { sum_{i in I(x*)} mu_i * nabla g_i(x*) | mu_i >= 0 }")

    # Build and print the equation with numerical values
    mu_vars = [f"mu_{i+1}" for i in active_indices]
    terms = [f"{m} * {g.T}" for m, g in zip(mu_vars, active_gradients)]
    equation = "s = " + " + ".join(terms)

    print("\nFor any vector s = (s_1, s_2, s_3)^T in T_F^°(x^*):")
    print(f"{equation}, where {', '.join(mu_vars)} >= 0.")
    
    print("\nWriting this component-wise:")
    print(f"s_1 = ({active_gradients[0][0]})*mu_1 + ({active_gradients[1][0]})*mu_2 = 2*mu_1 - 2*mu_2")
    print(f"s_2 = ({active_gradients[0][1]})*mu_1 + ({active_gradients[1][1]})*mu_2 = 0")
    print(f"s_3 = ({active_gradients[2][2]})*mu_3 = mu_3")

    print("\n--- Final Explicit Representation ---")
    print("Analyzing the components:")
    print("1. For s_1 = 2*(mu_1 - mu_2), since mu_1 >= 0 and mu_2 >= 0, s_1 can take any real value.")
    print("2. s_2 is always 0.")
    print("3. For s_3 = mu_3, since mu_3 >= 0, s_3 must be non-negative.")
    
    print("\nThus, the explicit representation of the normal cone is:")
    final_representation = "{ s = (s_1, s_2, s_3)^T in R^3 | s_2 = 0, s_3 >= 0 }"
    print(final_representation)

solve_normal_cone()