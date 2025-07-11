def solve():
    """
    This function demonstrates that for the torus T^2, the two actions of the
    fundamental group on the fiber of the universal cover are the same.
    """
    # An element of the fundamental group pi_1(T^2) ~= Z^2
    m = 2
    n = 3
    gamma_element = (m, n)

    # A point in the fiber p^{-1}(x_0) ~= Z^2
    k = 5
    l = 7
    fiber_point = (k, l)

    print(f"Let's consider the torus T^2.")
    print(f"The fundamental group pi_1(T^2) is isomorphic to Z^2.")
    print(f"Let's take a loop class represented by the element gamma = ({m}, {n}).")
    print(f"The universal cover is R^2, and the fiber over the basepoint is the integer lattice Z^2.")
    print(f"Let's take a point in the fiber, x_tilde = ({k}, {l}).")
    print("\nWe will now compute the result of the two actions on this point.")

    # --- Action 1: Action by restricting deck transformations ---
    print("\n--- Action 1: By Deck Transformations ---")
    print(f"The loop class gamma = ({m}, {n}) corresponds to a deck transformation T_gamma,")
    print(f"which is a translation by the vector ({m}, {n}).")
    print(f"The action on x_tilde = ({k}, {l}) is T_gamma(x_tilde).")
    
    result_1_k = k + m
    result_1_l = l + n
    print(f"Result = ({k} + {m}, {l} + {n}) = ({result_1_k}, {result_1_l})")

    # --- Action 2: Action by holonomy (path lifting) ---
    print("\n--- Action 2: By Holonomy (Path Lifting) ---")
    print(f"The loop for gamma = ({m}, {n}) is lifted to a path in the universal cover R^2,")
    print(f"starting at x_tilde = ({k}, {l}).")
    print(f"The lifted path is p(t) = ({m}*t + {k}, {n}*t + {l}) for t in [0, 1].")
    print(f"The action is the endpoint of this path (at t=1).")
    
    result_2_k = m * 1 + k
    result_2_l = n * 1 + l
    print(f"Result = ({m}*1 + {k}, {n}*1 + {l}) = ({result_2_k}, {result_2_l})")

    # --- Conclusion ---
    print("\n--- Conclusion ---")
    if (result_1_k, result_1_l) == (result_2_k, result_2_l):
        print("The results of the two actions are identical.")
    else:
        print("The results of the two actions are different.")

solve()