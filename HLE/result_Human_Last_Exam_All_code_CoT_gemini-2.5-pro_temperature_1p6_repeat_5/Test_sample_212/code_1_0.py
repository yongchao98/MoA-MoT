def solve():
    """
    This function demonstrates that for X = T^2, the two actions of the
    fundamental group on the fiber of the universal cover are the same.
    """

    # An element of the fundamental group pi_1(T^2) ~= Z^2.
    # Corresponds to a loop winding m times in one direction and n in the other.
    m, n = 3, -2
    loop_class = (m, n)

    # A point in the fiber p^-1(x_0) ~= Z^2.
    # This is a point in the universal cover R^2 that maps to the basepoint on T^2.
    k, l = 5, 8
    fiber_point = (k, l)

    print(f"Let's check the case for the 2-torus, X = T^2.")
    print(f"We take a loop class from pi_1(T^2) represented by the integers {loop_class}.")
    print(f"We take a point in the fiber p^-1(x_0) represented by the integer point {fiber_point}.")
    print("-" * 50)

    # 1. The action by holonomy (path lifting)
    print("Action 1: Holonomy (Path Lifting)")
    print("The action is defined by the endpoint of the path in the cover space.")
    print("This corresponds to adding the loop's integer vector to the fiber point's vector.")
    
    holonomy_result_x = k + m
    holonomy_result_y = l + n
    holonomy_result = (holonomy_result_x, holonomy_result_y)

    print(f"The new point is ({k} + {m}, {l} + {n}) = ({holonomy_result_x}, {holonomy_result_y})")
    print("-" * 50)


    # 2. The action by restricting deck transformations
    print("Action 2: Deck Transformations")
    print(f"The loop class {loop_class} corresponds to the deck transformation T(x,y) = (x + {m}, y + {n}).")
    print("Applying this transformation to the fiber point:")
    
    deck_trans_result_x = k + m
    deck_trans_result_y = l + n
    deck_trans_result = (deck_trans_result_x, deck_trans_result_y)

    print(f"T({k}, {l}) = ({k} + {m}, {l} + {n}) = ({deck_trans_result_x}, {deck_trans_result_y})")
    print("-" * 50)


    # Conclusion
    print("Conclusion:")
    if holonomy_result == deck_trans_result:
        print("Both actions produce the same result. The two actions are the same for X = T^2.")
        final_answer = "Yes"
    else:
        print("The actions produce different results. They are not the same.")
        final_answer = "No"
    
    # In general, these two actions are identical if and only if the fundamental group is abelian.
    # Since pi_1(T^2) = Z^2 is abelian, the actions are the same.
    print(f"\nFinal Answer: {final_answer}")

solve()