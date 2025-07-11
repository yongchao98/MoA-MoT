import sys

def calculate_actions_on_torus_fiber(loop_mn, fiber_point_ab):
    """
    This function simulates and compares two actions of the fundamental group
    of the 2-torus (T^2) on a point in the fiber of its universal cover.

    Args:
        loop_mn (tuple): A pair of integers (m, n) representing an element
                         of the fundamental group pi_1(T^2) ~= Z^2.
        fiber_point_ab (tuple): A pair of integers (a, b) representing a point
                                in the fiber, which is identified with Z^2.
    """

    m, n = loop_mn
    a, b = fiber_point_ab

    # We need a base point in the fiber to define the isomorphism between
    # the fundamental group and the group of deck transformations.
    # We choose the origin (0, 0).
    base_point_in_fiber = (0, 0)
    
    print("--- Setup ---")
    print(f"The chosen element of the fundamental group pi_1(T^2) corresponds to the integer vector (m, n) = {loop_mn}.")
    print(f"The chosen point in the fiber p^-1(x_0) is the integer vector (a, b) = {fiber_point_ab}.")
    print(f"The base point in the fiber is chosen to be {base_point_in_fiber}.")
    print("-" * 20)
    
    # 1. Action by Holonomy
    # This action is calculated by lifting the loop corresponding to (m, n)
    # to a path in the universal cover (R^2) starting at the fiber point (a, b).
    # The lift of the loop gamma_{m,n}(t) = (mt, nt) mod 1 starting at (a,b)
    # is the path tilde_gamma(t) = (mt + a, nt + b).
    # The action's result is the endpoint of this path at t=1.
    holonomy_result_x = a + m
    holonomy_result_y = b + n
    holonomy_result = (holonomy_result_x, holonomy_result_y)
    
    print("1. Calculating the action by holonomy:")
    print("   The result is the endpoint of the lift of the loop, which is (a+m, b+n).")
    print(f"   Calculation: ({a} + {m}, {b} + {n}) = {holonomy_result}")
    print("-" * 20)

    # 2. Action by Deck Transformations
    # First, find the deck transformation corresponding to the loop (m, n).
    # This is done by finding the endpoint of the lift of the loop starting
    # from the chosen base point in the fiber, (0, 0).
    # The lift from (0,0) ends at (m, n).
    # For the torus, the deck transformation that maps (0, 0) to (m, n)
    # is a translation by the vector (m, n).
    deck_translation_vector = (m, n)
    
    # Second, apply this deck transformation (translation) to the fiber point (a, b).
    deck_action_result_x = a + deck_translation_vector[0]
    deck_action_result_y = b + deck_translation_vector[1]
    deck_action_result = (deck_action_result_x, deck_action_result_y)

    print("2. Calculating the action by deck transformations:")
    print(f"   The deck transformation for loop {loop_mn} is a translation by the vector {deck_translation_vector}.")
    print("   Applying this translation to the fiber point (a,b) gives (a+m, b+n).")
    print(f"   Calculation: ({a} + {m}, {b} + {n}) = {deck_action_result}")
    print("-" * 20)

    # Compare the results
    are_same = (holonomy_result == deck_action_result)
    print("--- Conclusion ---")
    print(f"Result of holonomy action: {holonomy_result}")
    print(f"Result of deck transformation action: {deck_action_result}")
    print(f"Are the results of the two actions the same? {are_same}")

if __name__ == '__main__':
    # Example values for the loop class and the fiber point
    example_loop = (2, 3)
    example_fiber_point = (4, 5)
    
    # Run the simulation
    calculate_actions_on_torus_fiber(example_loop, example_fiber_point)