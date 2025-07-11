def solve_actions_on_torus():
    """
    This script demonstrates that for the 2-torus T^2, the action of the
    fundamental group on the fiber of the universal cover is the same whether
    defined by holonomy or by deck transformations.
    """

    # We represent points in the fiber p^{-1}(x_0) = Z^2 and elements of
    # the fundamental group pi_1(T^2, x_0) = Z^2 as tuples of integers.

    # Let's choose an arbitrary point in the fiber.
    # This corresponds to a point (m, n) in Z^2.
    p_tilde = (3, 5)

    # Let's choose an arbitrary element of the fundamental group.
    # This corresponds to a loop class represented by (a, b) in Z^2.
    g = (2, 7)

    # Unpack the tuples for clarity in the equations.
    m, n = p_tilde
    a, b = g

    print("--- Problem Setup ---")
    print(f"Space X = T^2 (the 2-torus)")
    print(f"Fundamental group pi_1(T^2) is represented by Z^2.")
    print(f"Fiber of the universal cover is also represented by Z^2.")
    print(f"Let's test with a point in the fiber p_tilde = ({m}, {n})")
    print(f"and an element of the fundamental group g = ({a}, {b}).")
    print("-" * 25)

    # --- Action 1: Holonomy ---
    # The result is the endpoint of the lift of the loop g starting from p_tilde.
    # For g=(a,b) and p_tilde=(m,n), this endpoint is (m+a, n+b).
    holonomy_result_x = m + a
    holonomy_result_y = n + b

    print("Action 1: By Holonomy")
    print("This action is calculated by lifting the loop g from the point p_tilde.")
    print(f"The resulting point's coordinates are given by the equation:")
    print(f"  ({m} + {a}, {n} + {b})")
    print(f"Result = ({holonomy_result_x}, {holonomy_result_y})")
    print("-" * 25)

    # --- Action 2: Deck Transformations ---
    # The element g=(a,b) corresponds to the deck transformation T_{a,b}.
    # Applying this transformation to p_tilde=(m,n) gives T_{a,b}(m,n) = (m+a, n+b).
    deck_result_x = m + a
    deck_result_y = n + b

    print("Action 2: By Deck Transformations")
    print("This action is calculated by applying the deck transformation corresponding to g.")
    print(f"The resulting point's coordinates are given by the equation:")
    print(f"  ({m} + {a}, {n} + {b})")
    print(f"Result = ({deck_result_x}, {deck_result_y})")
    print("-" * 25)

    # --- Conclusion ---
    print("Conclusion:")
    if (holonomy_result_x, holonomy_result_y) == (deck_result_x, deck_result_y):
        print("The results from both actions are identical.")
    else:
        print("The results from the actions are different.")

solve_actions_on_torus()