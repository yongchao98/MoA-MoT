def demonstrate_torus_actions():
    """
    Demonstrates that the two actions of the fundamental group of the torus T^2
    on the fiber of its universal cover are the same.
    """
    # An element of the fundamental group pi_1(T^2) is represented by a pair of integers.
    # This corresponds to a loop winding 'm' times along the first circle and 'n' times along the second.
    m, n = 3, 5
    loop_mn = (m, n)

    # A point in the fiber over the basepoint is also represented by a pair of integers.
    k, l = 1, 2
    fiber_point_kl = (k, l)

    print("We will check if the two actions of the fundamental group on the fiber are the same for X = T^2.")
    print("An element of the fundamental group pi_1(T^2) can be represented by (m, n).")
    print(f"A point in the fiber p^-1(x_0) can be represented by (k, l).")
    print(f"We choose an example: loop (m, n) = {loop_mn} acting on fiber point (k, l) = {fiber_point_kl}.\n")

    # --- Action 1: Holonomy ---
    print("Action 1: Holonomy around loops")
    print("---------------------------------")
    print("The action is defined by lifting the loop (m, n) to a path in the universal cover R^2, starting at the fiber point (k, l).")
    print("The result is the endpoint of this lifted path.")
    holonomy_result_k = k + m
    holonomy_result_l = l + n
    print("The final point's calculation is:")
    print(f"({k}, {l}) -> ({k} + {m}, {l} + {n}) = ({holonomy_result_k}, {holonomy_result_l})\n")


    # --- Action 2: Deck Transformations ---
    print("Action 2: Restricting deck transformations to the fiber")
    print("---------------------------------------------------------")
    print("First, we find the deck transformation corresponding to the loop (m, n).")
    print("This is the transformation that maps the basepoint of the fiber (0, 0) to the endpoint of the lift of the loop starting from (0, 0).")
    print(f"The endpoint of the lift from (0, 0) is ({m}, {n}). The deck transformation is phi(x, y) = (x + {m}, y + {n}).")
    print("\nNext, we apply this deck transformation to our fiber point (k, l).")
    deck_result_k = k + m
    deck_result_l = l + n
    print("The final point's calculation is:")
    print(f"phi({k}, {l}) = ({k} + {m}, {l} + {n}) = ({deck_result_k}, {deck_result_l})\n")

    # --- Conclusion ---
    print("Conclusion")
    print("----------")
    print("As the calculations show, both actions produce the exact same result.")
    print(f"Holonomy Action Result: ({holonomy_result_k}, {holonomy_result_l})")
    print(f"Deck Transformation Action Result: ({deck_result_k}, {deck_result_l})")
    print("\nThis holds for any integers m, n, k, l. For the torus T^2, the two actions are the same.")

demonstrate_torus_actions()