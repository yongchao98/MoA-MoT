def demonstrate_torus_actions():
    """
    Demonstrates that the two actions of the fundamental group on the fiber
    of the universal cover are the same for the 2-torus (T^2).
    """

    # Let's represent elements of the fundamental group pi_1(T^2) as integer pairs (a, b).
    # Let's represent points in the fiber p^-1(x_0) as integer pairs (m, n).
    
    # Example inputs:
    # An element of the fundamental group pi_1(T^2) ~= Z^2
    loop_element = (3, 2)
    # A point in the fiber p^-1(x_0) ~= Z^2
    fiber_point = (5, 7)

    a, b = loop_element
    m, n = fiber_point

    print("--- Setup for the 2-Torus (T^2) ---")
    print(f"Universal Cover (X_tilde): R^2")
    print(f"Fundamental Group (pi_1(T^2)): Z^2")
    print(f"Fiber (p^-1(x_0)): The integer lattice Z^2")
    print("-" * 35)
    print(f"Let's consider a loop representing the group element: [{a}, {b}]")
    print(f"And a point in the fiber: ({m}, {n})")
    print("-" * 35)

    # --- Action 1: Holonomy around loops ---
    # This action is defined by lifting the loop path from the fiber point
    # and finding the endpoint of the lifted path. For the torus, this
    # corresponds to vector addition.
    
    holonomy_result_x = m + a
    holonomy_result_y = n + b
    holonomy_result = (holonomy_result_x, holonomy_result_y)
    
    print("Action 1: By Holonomy (Path Lifting)")
    print("The action of the loop [{a}, {b}] on the fiber point ({m}, {n}) is calculated.")
    print(f"Result: ({m}, {n}) + ({a}, {b}) = ({m} + {a}, {n} + {b}) = {holonomy_result}")
    print()

    # --- Action 2: Restricting deck transformations ---
    # The deck transformation corresponding to the loop (a, b) is a
    # translation g(x, y) = (x+a, y+b). We apply this to the fiber point.
    
    deck_trans_result_x = m + a
    deck_trans_result_y = n + b
    deck_trans_result = (deck_trans_result_x, deck_trans_result_y)
    
    print("Action 2: By Restricting Deck Transformations")
    print(f"The deck transformation for loop [{a}, {b}] is g(x,y) = (x + {a}, y + {b}).")
    print(f"Applying this to the fiber point ({m}, {n}):")
    print(f"Result: g({m}, {n}) = ({m} + {a}, {n} + {b}) = {deck_trans_result}")
    print("-" * 35)

    # --- Conclusion ---
    if holonomy_result == deck_trans_result:
        print("Conclusion: The results are identical. The two actions are the same for T^2.")
    else:
        print("Conclusion: The results are different.")

if __name__ == '__main__':
    demonstrate_torus_actions()