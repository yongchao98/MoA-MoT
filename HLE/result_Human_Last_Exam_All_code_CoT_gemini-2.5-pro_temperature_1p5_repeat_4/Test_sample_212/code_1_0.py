def solve_torus_actions():
    """
    This script demonstrates that for the universal cover of the torus T^2,
    the action of the fundamental group on a fiber by holonomy is the same
    as the action by deck transformations.
    """

    # Let the space X be the torus T^2.
    # Its universal cover is the plane R^2, and the covering map is p(x,y) = (x mod 1, y mod 1).
    # The fundamental group pi_1(T^2) is isomorphic to Z^2.
    # An element of the fundamental group can be represented by a pair of integers (k, l),
    # corresponding to a loop that wraps k times longitudinally and l times latitudinally.
    group_element = (3, 2)  # An example element [gamma] from pi_1(T^2)

    # Let the basepoint in the torus be x_0 = (0, 0).
    # The fiber p^{-1}(x_0) is the set of points in R^2 that map to x_0.
    # This is the integer lattice Z^2.
    fiber_point = (5, 7)  # An example point ~y in the fiber p^{-1}(x_0)

    # We choose the basepoint in the universal cover ~x_0 to be (0, 0), which is in the fiber.
    
    print("Let's test the two actions for X = T^2.")
    print(f"Fundamental group element [gamma] corresponds to the vector (k, l) = {group_element}.")
    print(f"We will act on the fiber point ~y = (m, n) = {fiber_point}.")
    print("-" * 40)

    # --- Action 1: Holonomy around loops ---
    # This action is defined by lifting the loop gamma to a path in R^2 starting at ~y.
    # The result of the action is the endpoint of this lifted path.
    # Lifting a loop corresponding to (k, l) starting at (m, n) results in a path
    # from (m, n) to (m+k, n+l).
    k, l = group_element
    m, n = fiber_point
    
    result_holonomy_m = m + k
    result_holonomy_n = n + l
    result_holonomy = (result_holonomy_m, result_holonomy_n)
    
    print("1. Action by Holonomy ([gamma] . ~y):")
    print("   This is the endpoint of the lift of gamma starting at ~y.")
    print(f"   The result is ({m} + {k}, {n} + {l}) = ({result_holonomy_m}, {result_holonomy_n}).")
    print("-" * 40)


    # --- Action 2: Restricting deck transformations to the fiber ---
    # This action is defined by first finding the deck transformation g_gamma corresponding
    # to [gamma], and then applying it to the fiber point ~y.
    
    # The deck transformation g_gamma is the unique one that maps the basepoint ~x_0=(0,0)
    # to the endpoint of the lift of gamma starting at ~x_0.
    # The lift of the loop (k,l) starting at (0,0) ends at (k,l).
    # So, g_gamma must map (0,0) to (k,l). This means g_gamma is translation by (k,l).
    deck_translation_vector = group_element

    # Now, we apply this deck transformation (translation) to our fiber point ~y = (m, n).
    result_deck_m = m + deck_translation_vector[0]
    result_deck_n = n + deck_translation_vector[1]
    result_deck = (result_deck_m, result_deck_n)
    
    print("2. Action by Deck Transformation (g_gamma(~y)):")
    print(f"   The loop [gamma] corresponds to the deck transformation: translation by {deck_translation_vector}.")
    print(f"   Applying this translation to ~y = ({m}, {n}) gives:")
    print(f"   ({m} + {deck_translation_vector[0]}, {n} + {deck_translation_vector[1]}) = ({result_deck_m}, {result_deck_n}).")
    print("-" * 40)

    # --- Comparison ---
    print("Comparison:")
    print(f"Result from Holonomy Action:        {result_holonomy}")
    print(f"Result from Deck Transformation Action: {result_deck}")

    if result_holonomy == result_deck:
        print("\nThe two actions produce the same result.")
    else:
        print("\nThe two actions produce different results.")

solve_torus_actions()