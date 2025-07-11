import sys

def solve():
    """
    This script demonstrates that for the torus T^2, the action of the fundamental group
    on the fiber of the universal cover is the same whether defined by holonomy or
    by deck transformations.
    """

    # --- Setup ---
    # We model the mathematical objects using pairs of integers.
    # The universal cover of the torus T^2 is the plane R^2.
    # The fiber p^{-1}(x_0) is the set of integer lattice points, Z^2.
    # The fundamental group pi_1(T^2) is also isomorphic to Z^2.

    # Let's choose a sample point in the fiber p^{-1}(x_0).
    point_in_fiber = (5, 8)

    # Let's choose a sample element of the fundamental group pi_1(T^2).
    # This corresponds to a loop wrapping, say, 2 times in the first direction
    # and 3 times in the second direction.
    group_element = (2, 3)

    print(f"Demonstrating the two actions for X = T^2 (torus):")
    print(f"Let's consider a point in the fiber p^-1(x_0): {point_in_fiber}")
    print(f"And an element of the fundamental group pi_1(T^2): {group_element}")
    print("-" * 50)

    # --- Action 1: Holonomy around loops ---
    print("1. Action by Holonomy:")
    def action_by_holonomy(gamma, point):
        """
        Calculates the action of a fundamental group element on a fiber point via holonomy.
        This is done by lifting the loop corresponding to `gamma` to a path in the
        universal cover starting at `point`, and finding the endpoint.
        """
        m, n = point
        a, b = gamma
        # Lifting a loop corresponding to (a, b) from a starting point (m, n)
        # in the universal cover R^2 results in a path that ends at (m + a, n + b).
        result = (m + a, n + b)
        print(f"   Lifting the loop {gamma} starting from {point} ends at {result}.")
        print(f"   The final equation is: ({m}, {n}) + ({a}, {b}) = ({m + a}, {n + b})")
        return result

    result1 = action_by_holonomy(group_element, point_in_fiber)
    print("-" * 50)


    # --- Action 2: Restricting deck transformations ---
    print("2. Action by Deck Transformations:")
    def action_by_deck_transformation(gamma, point):
        """
        Calculates the action by finding the deck transformation corresponding
        to `gamma` and applying it to `point`.
        """
        m, n = point
        a, b = gamma

        # First, find the deck transformation g_gamma for the group element gamma = (a, b).
        # We identify g_gamma by its action on a basepoint in the fiber, say (0,0).
        # g_gamma((0,0)) must equal the endpoint of the lift of gamma from (0,0), which is (a,b).
        # Deck transformations for the torus are translations by integer vectors, g(x,y) = (x+c, y+d).
        # So g_gamma is the translation by the vector (a,b).
        # Now, we apply this transformation g_gamma(x,y) = (x+a, y+b) to our point (m, n).
        result = (m + a, n + b)
        print(f"   The deck transformation for loop {gamma} is a translation by the vector {gamma}.")
        print(f"   Applying this transformation to {point} gives {result}.")
        print(f"   The final equation is: g_({a},{b})(({m}, {n})) = ({m} + {a}, {n} + {b}) = ({m + a}, {n + b})")
        return result

    result2 = action_by_deck_transformation(group_element, point_in_fiber)
    print("-" * 50)


    # --- Conclusion ---
    print("Conclusion:")
    if result1 == result2:
        print("The results of the two actions are identical.")
        print("\nThis holds true because the fundamental group of the torus, pi_1(T^2) = Z^2, is abelian.")
        print("In general, these two actions are the same if and only if the fundamental group is abelian.")
    else:
        # This case should not be reached for the torus.
        print("The results are different.")

solve()