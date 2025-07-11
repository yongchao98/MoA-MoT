def solve_torus_action():
    """
    Demonstrates that for the 2-torus (T^2), the action of the fundamental
    group on a fiber by holonomy is the same as the action by restricting
    deck transformations.
    """
    # For T^2, both the fundamental group pi_1(T^2) and the fiber p^-1(x_0)
    # can be identified with the integer lattice Z^2.

    # Let's choose a sample element from the fundamental group.
    # This represents a loop going 2 times around one direction and 3 times
    # around the other.
    loop_g = (2, 3)

    # Let's choose a sample point in the fiber of the universal cover.
    fiber_point_x = (5, 7)

    # --- 1. Action by Holonomy (Path Lifting) ---
    def action_holonomy(loop, fiber_point):
        """
        Calculates the result of the holonomy action.
        Lifting the loop 'loop' starting from 'fiber_point' gives a new point
        which is the vector sum of the two.
        """
        k, l = loop
        m, n = fiber_point
        return (m + k, n + l)

    # --- 2. Action by Deck Transformations ---
    def action_deck(loop, fiber_point):
        """
        Calculates the result of the deck transformation action.
        The loop 'loop' corresponds to a deck transformation T(v) = v + loop.
        We apply this transformation to the 'fiber_point'.
        """
        k, l = loop
        m, n = fiber_point
        return (m + k, n + l)

    result_holonomy = action_holonomy(loop_g, fiber_point_x)
    result_deck = action_deck(loop_g, fiber_point_x)

    # --- Print the demonstration ---
    print(f"Let's choose an element from the fundamental group g = {loop_g}.")
    print(f"Let's choose a point in the fiber x_tilde = {fiber_point_x}.")
    print("-" * 60)

    print("1. Action by Holonomy (path lifting):")
    print("The new point is the endpoint of the path corresponding to g, lifted to start at x_tilde.")
    print("This corresponds to the equation: result = x_tilde + g")
    print(f"Calculation: {fiber_point_x} + {loop_g} = ({fiber_point_x[0]} + {loop_g[0]}, {fiber_point_x[1]} + {loop_g[1]}) = {result_holonomy}")
    print("-" * 60)

    print("2. Action by Deck Transformations:")
    print("The group element g corresponds to a deck transformation T_g(v) = v + g.")
    print("The action is applying T_g to x_tilde: result = T_g(x_tilde)")
    print(f"Calculation: T_{loop_g}({fiber_point_x}) = {fiber_point_x} + {loop_g} = ({fiber_point_x[0]} + {loop_g[0]}, {fiber_point_x[1]} + {loop_g[1]}) = {result_deck}")
    print("-" * 60)

    if result_holonomy == result_deck:
        print("Conclusion: The results are identical. The two actions are the same.")
    else:
        print("Conclusion: The results are different.")

solve_torus_action()