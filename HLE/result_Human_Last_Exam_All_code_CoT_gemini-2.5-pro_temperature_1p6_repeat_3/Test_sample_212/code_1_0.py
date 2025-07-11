def solve_torus_action():
    """
    Demonstrates the two actions of the fundamental group pi_1(T^2) on a fiber.

    For the 2-torus T^2, the universal cover is R^2, the fundamental group pi_1(T^2)
    is isomorphic to Z^2, and the fiber p^-1(x_0) is also identifiable with Z^2.

    This script shows that the action by holonomy and the action by deck transformations
    are equivalent for T^2.
    """

    # We represent an element of the fundamental group pi_1(T^2) as an integer tuple (m, n).
    # This corresponds to a loop winding m times in one direction and n times in the other.
    gamma_loop_class = (3, 5)

    # We represent a point in the fiber p^-1(x_0) as an integer tuple (k, l).
    fiber_point = (2, 7)

    m, n = gamma_loop_class
    k, l = fiber_point

    print(f"Let's consider the action of the loop class [gamma] = {gamma_loop_class} on the fiber point x_tilde = {fiber_point}.\n")

    # ----- Action 1: Holonomy around loops -----
    # The action is defined by the endpoint of the lift of the loop 'gamma' starting from 'fiber_point'.
    # A path for the loop is gamma(t) = (m*t, n*t) mod 1.
    # The lift of this path starting at (k, l) is tilde_gamma(t) = (m*t + k, n*t + l).
    # The endpoint at t=1 is (m*1 + k, n*1 + l).
    holonomy_result_x = m + k
    holonomy_result_y = n + l
    holonomy_result = (holonomy_result_x, holonomy_result_y)

    print("1. Action by Holonomy (Path Lifting):")
    print(f"   The lifted path starts at {fiber_point} and ends at the new point.")
    print(f"   The final equation for the new point is: ([gamma] acted on x_tilde) = (m+k, n+l)")
    print(f"   Substituting the numbers: ({m} + {k}, {n} + {l}) = {holonomy_result}\n")

    # ----- Action 2: Restricting deck transformations to the fiber -----
    # The loop class [gamma] = (m, n) corresponds to the deck transformation g(x, y) = (x+m, y+n).
    # This transformation is applied to the fiber point.
    deck_result_x = k + m
    deck_result_y = l + n
    deck_result = (deck_result_x, deck_result_y)

    print("2. Action by Deck Transformations:")
    print(f"   The loop class {gamma_loop_class} corresponds to the deck transformation g(x, y) = (x+{m}, y+{n}).")
    print(f"   Applying this transformation to the fiber point {fiber_point}:")
    print(f"   The final equation is: g({k}, {l}) = ({k}+{m}, {l}+{n})")
    print(f"   Substituting the numbers: ({k} + {m}, {l} + {n}) = {deck_result}\n")
    
    # ----- Conclusion -----
    print("="*40)
    print("Conclusion:")
    if holonomy_result == deck_result:
        print("The results of both actions are identical.")
        print("This holds true in general for the torus because its fundamental group, Z^2, is abelian.")
    else:
        # This case will not be reached for the torus.
        print("The results of the actions are different.")

solve_torus_action()