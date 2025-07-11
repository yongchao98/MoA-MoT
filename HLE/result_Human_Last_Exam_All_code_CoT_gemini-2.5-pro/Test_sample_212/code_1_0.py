def solve_torus_actions():
    """
    This script demonstrates the two actions of the fundamental group pi_1(T^2)
    on a fiber of its universal cover and shows they are identical.
    """
    # We will use an example to illustrate the actions.
    # Let the element of the fundamental group pi_1(T^2) ~= Z^2 be:
    gamma = (2, 3)
    # Let the point in the fiber p^{-1}(x_0) ~= Z^2 be:
    fiber_point = (5, 7)

    a, b = gamma
    m, n = fiber_point

    print("We are examining two actions of the fundamental group pi_1(T^2) on a fiber of its universal cover.")
    print(f"Let's use the group element gamma = {gamma} and the fiber point = {fiber_point}.\n")

    # --- Action 1: Holonomy (Path Lifting) ---
    print("--- 1. Action by Holonomy (Path Lifting) ---")
    print("This action is defined by lifting the loop corresponding to gamma to a path starting at the fiber point.")
    print(f"The lift of the loop for gamma={gamma} starting at {fiber_point} is the path p_tilde(t) = ({a}*t + {m}, {b}*t + {n}).")
    print("The result of the action is the endpoint of this path at t=1.")
    
    result_holonomy_x = m + a
    result_holonomy_y = n + b
    
    print("The final equation for the holonomy action is:")
    print(f"gamma . fiber_point = ({a}, {b}) . ({m}, {n}) = ({m} + {a}, {n} + {b}) = ({result_holonomy_x}, {result_holonomy_y})")

    # --- Action 2: Deck Transformations ---
    print("\n--- 2. Action by Deck Transformations ---")
    print("This action identifies the group element with a deck transformation and applies it.")
    print(f"The deck transformation f_gamma for gamma={gamma} is the translation by the vector {gamma}.")
    print(f"We apply this translation to the fiber point {fiber_point}.")

    result_deck_x = m + a
    result_deck_y = n + b

    print("The final equation for the deck transformation action is:")
    print(f"f_gamma(fiber_point) = T_({a},{b})({m},{n}) = ({m} + {a}, {n} + {b}) = ({result_deck_x}, {result_deck_y})")

    # --- Conclusion ---
    print("\n--- Conclusion ---")
    if (result_holonomy_x, result_holonomy_y) == (result_deck_x, result_deck_y):
        print("The results of the two actions are identical.")
        print("For X = T^2, the two actions are indeed the same.")
    else:
        print("The results of the two actions are different.")

solve_torus_actions()