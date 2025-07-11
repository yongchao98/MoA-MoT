def demonstrate_torus_actions():
    """
    This function demonstrates that for the torus T^2, the action of the fundamental
    group on a fiber of the universal cover is the same, whether defined by
    holonomy or by deck transformations.

    Let X = T^2, the 2-torus.
    - The universal cover is tilde(X) = R^2.
    - The covering map is p(x, y) = (x mod 1, y mod 1).
    - Let the base point be x_0 = p(0,0) on T^2.
    - The fiber over x_0 is p^{-1}(x_0), which is the integer lattice Z^2 in R^2.
    - The fundamental group pi_1(T^2, x_0) is isomorphic to Z^2.
    """
    # Define an example element from the fundamental group pi_1(T^2) ~= Z^2
    group_element = (2, 7)
    m, n = group_element

    # Define an example point from the fiber p^{-1}(x_0) ~= Z^2
    fiber_point = (3, 5)
    k, l = fiber_point

    print(f"We will check the two actions of the fundamental group on the fiber for X = T^2.")
    print(f"Let's take the group element [{m}, {n}] from pi_1(T^2) ~= Z^2.")
    print(f"And the fiber point ({k}, {l}) from the fiber p^-1(x_0) ~= Z^2.\n")

    # --- Action 1: By holonomy (path lifting) ---
    print("--- Action 1: Holonomy ---")
    # A loop representing (m, n) in pi_1 is one whose lift from (0,0) ends at (m,n).
    # Lifting this same loop from a different point (k, l) in the fiber
    # results in a path from (k, l) to (k+m, l+n). The action is the endpoint.
    result_k_holo = k + m
    result_l_holo = l + n

    print(f"The action of [{m}, {n}] on ({k}, {l}) is the endpoint of the lifted path.")
    print(f"The resulting point is given by the equation: ({k}) + ({m}), ({l}) + ({n})")
    print(f"Resulting point: ({result_k_holo}, {result_l_holo})")
    print(f"Equation evaluated: {k} + {m} = {result_k_holo}, {l} + {n} = {result_l_holo}\n")


    # --- Action 2: By restricting deck transformations ---
    print("--- Action 2: Deck Transformations ---")
    # The group element (m, n) corresponds to the deck transformation g(x,y) = (x+m, y+n).
    # Applying this transformation to the fiber point (k, l):
    result_k_deck = k + m
    result_l_deck = l + n

    print(f"The action of [{m}, {n}] on ({k}, {l}) is applying the deck transformation g(x,y)=(x+{m},y+{n}).")
    print(f"The resulting point is given by the equation: g({k}, {l}) = ({k} + {m}, {l} + {n})")
    print(f"Resulting point: ({result_k_deck}, {result_l_deck})")
    print(f"Equation evaluated: {k} + {m} = {result_k_deck}, {l} + {n} = {result_l_deck}\n")

    # --- Conclusion ---
    print("--- Conclusion ---")
    if (result_k_holo, result_l_holo) == (result_k_deck, result_l_deck):
        print("For X = T^2, the results of the two actions are identical.")
        print("This is because the fundamental group of the torus, Z^2, is abelian, which makes the two definitions coincide.")
    else:
        # This case will not be reached for the torus.
        print("The results of the two actions are different.")

demonstrate_torus_actions()