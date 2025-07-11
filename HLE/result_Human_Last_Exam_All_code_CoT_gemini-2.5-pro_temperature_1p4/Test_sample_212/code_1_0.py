def demonstrate_torus_actions(loop_class, fiber_point):
    """
    Demonstrates the two actions of the fundamental group of a torus
    on the fiber of its universal cover.

    Args:
        loop_class (tuple): A pair of integers (k, l) representing the
                            loop class [alpha^k * beta^l] in pi_1(T^2).
        fiber_point (tuple): A pair of integers (m, n) representing a
                             point in the fiber Z^2.
    """
    k, l = loop_class
    m, n = fiber_point

    print(f"Let's consider the torus T^2.")
    print(f"The fundamental group pi_1(T^2) is isomorphic to Z^2.")
    print(f"An element of this group is represented by integers (k, l), for this example we use {loop_class}.")
    print(f"The universal cover is R^2, and the fiber over the basepoint is the integer lattice Z^2.")
    print(f"We will see how the loop acts on a point in the fiber, (m, n) = {fiber_point}.\n")

    # Action 1: Holonomy (path lifting)
    holo_result_m = m + k
    holo_result_n = n + l
    holo_result = (holo_result_m, holo_result_n)

    print("--- Action 1: Holonomy (Path Lifting) ---")
    print(f"To find the action of the loop class ({k}, {l}) on the fiber point ({m}, {n}), we lift the loop to a path in R^2 starting at ({m}, {n}).")
    print(f"This is equivalent to adding the vector representing the loop class to the vector of the fiber point.")
    print(f"The endpoint of the lifted path is calculated as:")
    print(f"({m} + {k}, {n} + {l}) = ({holo_result_m}, {holo_result_n})")
    print(f"Result of Action 1: {holo_result}\n")

    # Action 2: Deck Transformations
    deck_result_m = m + k
    deck_result_n = n + l
    deck_result = (deck_result_m, deck_result_n)

    print("--- Action 2: Deck Transformations ---")
    print(f"The deck transformation corresponding to the loop class ({k}, {l}) is the translation map T(x, y) = (x + {k}, y + {l}).")
    print(f"We apply this transformation to the fiber point ({m}, {n}):")
    print(f"T({m}, {n}) = ({m} + {k}, {n} + {l}) = ({deck_result_m}, {deck_result_n})")
    print(f"Result of Action 2: {deck_result}\n")

    # Conclusion
    are_same = (holo_result == deck_result)
    print("--- Conclusion ---")
    if are_same:
        print(f"The results {holo_result} and {deck_result} are identical.")
        print("Therefore, for the torus T^2, the two actions are the same.")
    else:
        # This case is not expected for the torus.
        print(f"The results {holo_result} and {deck_result} are different.")
        print("The two actions are not the same.")

# Let's run the demonstration with an example.
# We choose a loop class corresponding to going twice around the first circle
# and three times around the second circle of the torus.
example_loop = (2, 3)
# We choose an arbitrary point in the fiber.
example_point = (5, 7)
demonstrate_torus_actions(example_loop, example_point)