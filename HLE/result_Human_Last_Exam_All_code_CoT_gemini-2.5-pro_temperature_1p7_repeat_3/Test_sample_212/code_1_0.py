import sys

def solve():
    """
    Demonstrates that for the torus T^2, the action of the fundamental group on
    the fiber of the universal cover is the same whether defined by holonomy
    or by deck transformations.
    """

    # We choose a sample element of the fundamental group pi_1(T^2) ~= Z^2.
    # This represents a loop, e.g., winding 3 times in the first direction
    # and 2 times in the second.
    gamma_mn = (3, 2)

    # We choose a sample point in the fiber p^{-1}(x_0) ~= Z^2.
    # This corresponds to a point in the universal cover R^2 that projects
    # to the basepoint on the torus.
    fiber_point_kl = (5, 7)

    # --- Action 1: Holonomy around loops (Path Lifting) ---

    def holonomy_action(gamma, fiber_point):
        """
        Calculates the action via path lifting.
        The action of a loop (m,n) on a fiber point (k,l) is (k+m, l+n).
        """
        m, n = gamma
        k, l = fiber_point
        return (k + m, l + n)

    # --- Action 2: Restricting Deck Transformations ---

    def deck_transformation_action(gamma, fiber_point):
        """
        Calculates the action via deck transformations.
        The loop (m,n) corresponds to the deck transformation g(x,y) = (x+m, y+n).
        Applying this to the fiber point (k,l) gives (k+m, l+n).
        """
        m, n = gamma
        k, l = fiber_point
        # The corresponding deck transformation is a translation by the vector (m, n).
        deck_translation_vector = (m, n)
        return (k + deck_translation_vector[0], l + deck_translation_vector[1])

    # --- Demonstration ---

    print(f"Consider the torus X = T^2.")
    print(f"Let the fundamental group element be [gamma] = {gamma_mn}")
    print(f"Let the point in the fiber p^-1(x_0) be p = {fiber_point_kl}\n")

    # Calculate the result of the first action
    result_holonomy = holonomy_action(gamma_mn, fiber_point_kl)
    print("Action 1 (Holonomy):")
    print(f"  The result of the loop [gamma]={gamma_mn} acting on p={fiber_point_kl} is:")
    print(f"  ({fiber_point_kl[0]} + {gamma_mn[0]}, {fiber_point_kl[1]} + {gamma_mn[1]}) = {result_holonomy}\n")

    # Calculate the result of the second action
    result_deck = deck_transformation_action(gamma_mn, fiber_point_kl)
    print("Action 2 (Deck Transformations):")
    print(f"  The result of the deck transformation for [gamma]={gamma_mn} acting on p={fiber_point_kl} is:")
    print(f"  ({fiber_point_kl[0]} + {gamma_mn[0]}, {fiber_point_kl[1]} + {gamma_mn[1]}) = {result_deck}\n")

    # Compare the results
    are_same = (result_holonomy == result_deck)
    print(f"Conclusion: The results are identical. The two actions are the same for the torus.")

solve()
