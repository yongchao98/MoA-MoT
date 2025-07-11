import sys

def solve_actions_on_torus():
    """
    This script investigates two different actions of the fundamental group pi_1(T^2)
    on the fiber of its universal cover and determines if they are the same.
    """

    # Let's represent elements from the fundamental group pi_1(T^2) ~= Z^2 and
    # the fiber p^-1(x_0) ~= Z^2 using tuples of integers.
    # We choose an arbitrary element from the fundamental group.
    m, n = 3, 5
    gamma = (m, n)

    # We choose an arbitrary point in the fiber.
    k, l = 2, 7
    fiber_point = (k, l)

    print("--- Analyzing two group actions for the Torus T^2 ---")
    print(f"Let the fundamental group element [gamma] be represented by the integer pair ({m}, {n}).")
    print(f"Let the point in the fiber x_tilde be represented by the integer pair ({k}, {l}).")
    print("-" * 55)

    # 1. Action by restricting deck transformations
    print("Action 1: By restricting deck transformations to the fiber")

    # For the universal cover R^2 -> T^2, the element (m, n) of the fundamental group
    # corresponds to the deck transformation that is a translation by the vector (m, n).
    deck_translation_vector = gamma
    
    # We apply this translation to the fiber point.
    result1_x = fiber_point[0] + deck_translation_vector[0]
    result1_y = fiber_point[1] + deck_translation_vector[1]
    result1 = (result1_x, result1_y)
    
    print(f"The deck transformation for [gamma]=({m},{n}) is a translation by the vector ({m},{n}).")
    print(f"Applying this translation to the fiber point ({k},{l}):")
    # Output the final equation with each number
    print(f"  ({fiber_point[0]}, {fiber_point[1]}) + ({deck_translation_vector[0]}, {deck_translation_vector[1]}) = ({result1[0]}, {result1[1]})")
    print(f"Result from Action 1 is the point: {result1}")
    print("-" * 55)

    # 2. Action by holonomy around loops
    print("Action 2: By holonomy around loops in X")
    
    # The action of [gamma]=(m,n) on a fiber point (k,l) is the endpoint of the
    # lift of the loop gamma, where the lift starts at (k,l).
    # The lift of the path (t*m, t*n) starting at (k,l) is (t*m + k, t*n + l).
    # The endpoint is at t=1.
    result2_x = m + k
    result2_y = n + l
    result2 = (result2_x, result2_y)

    print(f"We lift the loop [gamma]=({m},{n}) starting from the fiber point ({k},{l}).")
    print(f"The endpoint of this lifted path is calculated as:")
    # Output the final equation with each number
    print(f"  ({k}, {l}) + ({m}, {n}) = ({result2[0]}, {result2[1]})")
    print(f"Result from Action 2 is the point: {result2}")
    print("-" * 55)

    # 3. Compare the results
    print("Comparison:")
    print(f"Result 1 (Deck Transformations): {result1}")
    print(f"Result 2 (Holonomy):            {result2}")
    
    if result1 == result2:
        print("\nThe two actions produce the same result.")
        final_answer = "Yes"
    else:
        # This case will not be reached for T^2.
        print("\nThe two actions produce different results.")
        final_answer = "No"

    # Final conclusion based on the demonstration. The mathematical reasoning is that
    # the two actions are the same if and only if the group of deck transformations
    # (which is isomorphic to the fundamental group) is abelian. For T^2, the group
    # is Z^2, which is abelian.
    sys.stdout.write("<<<" + final_answer + ">>>")

solve_actions_on_torus()