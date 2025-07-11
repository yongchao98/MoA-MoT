def solve_actions_on_torus():
    """
    Demonstrates that for the torus T^2, the action of the fundamental group on
    the fiber of the universal cover is the same whether defined by holonomy
    or by deck transformations.
    """
    # Let X = T^2. We identify its fundamental group pi_1(X) with Z^2.
    # We identify the fiber p^-1(x_0) over the basepoint with Z^2.

    # 1. Define a sample point in the fiber Z^2
    k = 5
    l = 7
    fiber_point = (k, l)

    # 2. Define a sample element from the fundamental group Z^2
    m = 2
    n = 3
    group_element = (m, n)

    print("--- Verifying actions of pi_1(T^2) on the fiber p^-1(x_0) ---")
    print(f"Identifying both pi_1(T^2) and the fiber with Z^2.")
    print(f"Sample fiber point (k, l): {fiber_point}")
    print(f"Sample group element (m, n): {group_element}\n")

    # 3. Calculate the result of Action 1 (Holonomy)
    # The result is the endpoint of the lift of the loop (m,n) starting at (k,l).
    # Lifted path: t -> (k + m*t, l + n*t)
    # Endpoint at t=1: (k + m, l + n)
    holo_res_k = k + m
    holo_res_l = l + n
    holonomy_result = (holo_res_k, holo_res_l)

    print("--- Action 1: By Holonomy ---")
    print("This action is calculated by lifting the loop and finding the endpoint.")
    print(f"Lifting the loop corresponding to {group_element} starting at {fiber_point} gives a path whose endpoint is:")
    print(f"({k} + {m}, {l} + {n}) = {holonomy_result}")
    print("The final equation is:")
    print(f"{k} + {m} = {holo_res_k}")
    print(f"{l} + {n} = {holo_res_l}\n")

    # 4. Calculate the result of Action 2 (Deck Transformations)
    # The group element (m,n) corresponds to the deck transformation g(x,y) = (x+m, y+n).
    # We apply this transformation to the fiber point (k,l).
    deck_res_k = k + m
    deck_res_l = l + n
    deck_result = (deck_res_k, deck_res_l)

    print("--- Action 2: By Deck Transformations ---")
    print("This action is calculated by applying the corresponding deck transformation.")
    print(f"The deck transformation for {group_element} is g(x,y) = (x+{m}, y+{n}).")
    print(f"Applying this to {fiber_point}: g({k}, {l}) = ({k}+{m}, {l}+{n}) = {deck_result}")
    print("The final equation is:")
    print(f"{k} + {m} = {deck_res_k}")
    print(f"{l} + {n} = {deck_res_l}\n")
    
    # 5. Compare the results
    are_same = (holonomy_result == deck_result)
    print("--- Comparison ---")
    print(f"Result from Holonomy Action:      {holonomy_result}")
    print(f"Result from Deck Transform Action: {deck_result}")
    print(f"Are the two actions the same for X=T^2? {are_same}")


if __name__ == '__main__':
    solve_actions_on_torus()