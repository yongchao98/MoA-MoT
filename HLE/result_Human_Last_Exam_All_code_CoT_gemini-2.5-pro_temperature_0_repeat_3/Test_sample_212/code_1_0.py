def demonstrate_torus_actions():
    """
    Demonstrates that for the torus T^2, the action by holonomy and the action
    by deck transformations on the fiber of the universal cover are the same.
    """
    # Let's represent an element of the fundamental group pi_1(T^2) ~= Z^2 as a tuple (m, n).
    # This corresponds to a loop winding m times longitudinally and n times meridionally.
    gamma = (3, 5)

    # Let's represent a point in the fiber p^-1(x_0) ~= Z^2 as a tuple (k, l).
    # These are the points in the universal cover R^2 that project to the basepoint on the torus.
    fiber_point = (2, 7)

    # --- Action 1: Holonomy around loops ---
    # The action is defined by the endpoint of the lift of the loop `gamma` starting from `fiber_point`.
    # Lifting a path corresponding to (m, n) from a start point (k, l) in R^2
    # results in a path ending at (k+m, l+n).
    def action_holonomy(gamma_mn, point_kl):
        m, n = gamma_mn
        k, l = point_kl
        return (k + m, l + n)

    # --- Action 2: Restricting deck transformations to the fiber ---
    # The deck transformation corresponding to gamma = (m, n) is phi(x, y) = (x+m, y+n).
    # The action is to apply this transformation to the `fiber_point`.
    def action_deck_transformation(gamma_mn, point_kl):
        m, n = gamma_mn
        k, l = point_kl
        # Apply the deck transformation phi(x,y) = (x+m, y+n) to the point (k,l)
        return (k + m, l + n)

    # --- Verification ---
    result1 = action_holonomy(gamma, fiber_point)
    result2 = action_deck_transformation(gamma, fiber_point)

    print("Let X = T^2 (the 2-torus).")
    print("The universal cover is R^2 and the fundamental group pi_1(T^2) is Z^2.")
    print("The fiber p^-1(x_0) over the basepoint is the integer lattice Z^2.")
    print("-" * 50)
    print(f"Let's consider a loop class gamma = {gamma} in pi_1(T^2)")
    print(f"and a point p = {fiber_point} in the fiber.")
    print("-" * 50)

    print("Action 1 (Holonomy):")
    print(f"The action of gamma on p is the endpoint of the lifted loop.")
    print(f"Result: {gamma} * {fiber_point} = ({fiber_point[0]} + {gamma[0]}, {fiber_point[1]} + {gamma[1]}) = {result1}")
    print()

    print("Action 2 (Deck Transformations):")
    print(f"The deck transformation for gamma = {gamma} is phi(x, y) = (x + {gamma[0]}, y + {gamma[1]}).")
    print(f"Applying this to p: phi{fiber_point} = ({fiber_point[0]} + {gamma[0]}, {fiber_point[1]} + {gamma[1]}) = {result2}")
    print()

    print("-" * 50)
    if result1 == result2:
        print("Conclusion: The results of the two actions are identical.")
    else:
        print("Conclusion: The results of the two actions are different.")

demonstrate_torus_actions()