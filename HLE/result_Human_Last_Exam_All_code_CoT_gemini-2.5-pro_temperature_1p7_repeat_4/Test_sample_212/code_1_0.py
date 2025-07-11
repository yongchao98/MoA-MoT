def demonstrate_torus_actions(group_element, fiber_point):
    """
    Demonstrates that the two actions of the fundamental group on the fiber
    of the universal cover are the same for the 2-torus.

    Args:
        group_element (tuple): A pair of integers (m, n) representing an element
                               of the fundamental group pi_1(T^2) = Z^2.
        fiber_point (tuple): A pair of integers (a, b) representing a point
                             in the fiber Z^2.
    """
    m, n = group_element
    a, b = fiber_point

    print("--- Demonstrating Actions on the Torus Fiber ---")
    print(f"Chosen fundamental group element: g = ({m}, {n})")
    print(f"Chosen point in the fiber p^{-1}(x_0): p = ({a}, {b})")
    print("--------------------------------------------------")

    # Action 1: Holonomy (Path Lifting)
    print("Action 1: By holonomy (path lifting)")
    # The result is the endpoint of the lift of the loop for (m, n), starting at (a, b).
    # The lifted path is L(t) = (a + m*t, b + n*t).
    # The endpoint L(1) is (a + m, b + n).
    holonomy_result_a = a + m
    holonomy_result_b = b + n
    print(f"The result is the endpoint of the lifted path:")
    print(f"L(1) = ({a} + {m}, {b} + {n}) = ({holonomy_result_a}, {holonomy_result_b})")
    print("--------------------------------------------------")

    # Action 2: Deck Transformations
    print("Action 2: By restricting deck transformations")
    # The deck transformation T_g for g=(m,n) is translation by the vector (m,n).
    # T_g(x, y) = (x + m, y + n).
    # We apply this transformation to the fiber point p=(a,b).
    deck_result_a = a + m
    deck_result_b = b + n
    print(f"The deck transformation for g=({m},{n}) is translation by ({m},{n}).")
    print(f"Applying this transformation to p=({a},{b}):")
    print(f"T_g(p) = ({a} + {m}, {b} + {n}) = ({deck_result_a}, {deck_result_b})")
    print("--------------------------------------------------")

    # Comparison
    print("Comparison:")
    if (holonomy_result_a, holonomy_result_b) == (deck_result_a, deck_result_b):
        print("The results are identical. The two actions are the same for the torus.")
    else:
        print("The results are different.")

# --- Execute the demonstration with an example ---
# Let's choose a non-trivial group element and fiber point.
pi_1_element = (2, -3)
fiber_point_element = (5, 1)
demonstrate_torus_actions(pi_1_element, fiber_point_element)