def solve_torus_actions():
    """
    Demonstrates the two actions of the fundamental group pi_1(T^2) on a fiber
    of its universal cover.
    """

    # We are considering the case where X is the 2-torus, T^2.
    # The fundamental group, pi_1(T^2), is isomorphic to Z^2.
    # We represent an element of the group, a loop class gamma, as an integer tuple (m, n).
    m = 2
    n = 3
    gamma = (m, n)

    # The universal cover X_tilde is R^2.
    # A point in the fiber p^-1(x_0) is a point on the integer lattice Z^2.
    # We represent such a point, x_tilde, as an integer tuple (a, b).
    a = 5
    b = 7
    x_tilde = (a, b)
    
    # We choose a base point in the fiber to establish the isomorphism between
    # pi_1 and the group of deck transformations. A standard choice is (0, 0).
    base_point_in_fiber = (0, 0)

    print(f"Let's consider the torus X = T^2.")
    print(f"The fundamental group element is gamma = {gamma}.")
    print(f"The point in the fiber is x_tilde = {x_tilde}.\n")

    # --- Action 1: Holonomy around loops (Monodromy action) ---
    # This action is computed by lifting the loop gamma to a path in the universal cover
    # starting at the fiber point x_tilde. The result is the endpoint of this path.
    # A lift of the loop corresponding to (m, n) starting at (a, b) ends at (a+m, b+n).
    
    holonomy_result_x = x_tilde[0] + gamma[0]
    holonomy_result_y = x_tilde[1] + gamma[1]
    holonomy_result = (holonomy_result_x, holonomy_result_y)

    print("Action 1: Holonomy around loops")
    print("This action lifts the loop gamma to a path starting at x_tilde.")
    print(f"The endpoint of this lifted path is the result of the action.")
    print(f"Calculation: ({a} + {m}, {b} + {n}) = {holonomy_result}\n")


    # --- Action 2: Restricting deck transformations to the fiber ---
    # First, we identify the deck transformation g_gamma that corresponds to gamma.
    # For T^2, deck transformations are translations by integer vectors.
    # The loop gamma=(m,n) corresponds to the translation by the vector (m,n).
    # So, g_gamma(x, y) = (x + m, y + n).
    # We then apply this transformation to our fiber point x_tilde.

    deck_result_x = x_tilde[0] + gamma[0]
    deck_result_y = x_tilde[1] + gamma[1]
    deck_result = (deck_result_x, deck_result_y)
    
    print("Action 2: Restricting deck transformations")
    print("This action finds the deck transformation g_gamma corresponding to gamma and applies it to x_tilde.")
    print(f"The transformation g_gamma is a translation by the vector {gamma}.")
    print(f"Applying this to x_tilde gives the result.")
    print(f"Calculation: ({a} + {m}, {b} + {n}) = {deck_result}\n")
    

    # --- Conclusion ---
    print("Conclusion:")
    if holonomy_result == deck_result:
        print("The results of the two actions are identical.")
        print("Therefore, for X = T^2, the two actions are the same.")
    else:
        # This case should not be reached for the torus.
        print("The results of the two actions are different.")

solve_torus_actions()