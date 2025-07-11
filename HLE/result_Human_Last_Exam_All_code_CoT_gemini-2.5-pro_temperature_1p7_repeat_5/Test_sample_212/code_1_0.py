def solve_actions_on_torus():
    """
    Demonstrates and compares the two actions of the fundamental group
    pi_1(T^2) on a fiber of its universal cover.
    """

    # === Setup ===
    # Let X be the torus T^2.
    # The fundamental group pi_1(T^2) is isomorphic to Z^2.
    # A loop gamma is represented by a pair of integers (m, n),
    # corresponding to the loop alpha^m * beta^n.
    gamma_m = 2
    gamma_n = 3
    gamma = (gamma_m, gamma_n)

    # The universal cover is R^2. The basepoint in X is x0 = [(0,0)].
    # The fiber p^-1(x0) is the integer lattice Z^2.
    # Let's choose a point y in the fiber.
    y_x = 5
    y_y = 7
    y = (y_x, y_y)

    print(f"Let X = T^2. We consider the action of a loop gamma = alpha^{gamma_m} * beta^{gamma_n} from pi_1(T^2)")
    print(f"on a point y = ({y_x}, {y_y}) in the fiber p^-1(x0).")
    print("-" * 50)

    # === Action 1: Holonomy ===
    print("Action 1: By holonomy (path lifting)")
    # We lift the path gamma starting from y. The action of a group element
    # [alpha^m beta^n] is calculated by lifting alpha^m and then beta^n.
    # Lifting alpha^m from y corresponds to adding (m, 0).
    intermediate_point_x = y_x + gamma_m
    intermediate_point_y = y_y
    intermediate_point = (intermediate_point_x, intermediate_point_y)
    print(f"Lifting alpha^{gamma_m} from y = {y} moves the point to ({y_x} + {gamma_m}, {y_y}) = {intermediate_point}.")

    # From this new point, lifting beta^n corresponds to adding (0, n).
    final_point1_x = intermediate_point_x
    final_point1_y = intermediate_point_y + gamma_n
    final_point1 = (final_point1_x, final_point1_y)
    print(f"Lifting beta^{gamma_n} from {intermediate_point} moves the point to ({intermediate_point_x}, {intermediate_point_y} + {gamma_n}) = {final_point1}.")
    print(f"Final Result of Action 1: {final_point1}")
    print("-" * 50)


    # === Action 2: Deck Transformations ===
    print("Action 2: By restricting deck transformations")
    # First, find the deck transformation g_gamma corresponding to gamma = (m,n).
    # The isomorphism is based at x_tilde_0 = (0,0). g_gamma maps (0,0) to the
    # endpoint of the lift of gamma from (0,0), which is (m,n).
    # For T^2, this transformation is a simple translation.
    g_gamma_translation_x = gamma_m
    g_gamma_translation_y = gamma_n
    print(f"The deck transformation for gamma = ({gamma_m}, {gamma_n}) is a translation by the vector ({g_gamma_translation_x}, {g_gamma_translation_y}).")

    # Now, apply this transformation to the point y.
    final_point2_x = y_x + g_gamma_translation_x
    final_point2_y = y_y + g_gamma_translation_y
    final_point2 = (final_point2_x, final_point2_y)
    print(f"Applying this translation to y = {y}: ({y_x} + {g_gamma_translation_x}, {y_y} + {g_gamma_translation_y}) = {final_point2}.")
    print(f"Final Result of Action 2: {final_point2}")
    print("-" * 50)

    # === Comparison ===
    print("Comparing the results:")
    print(f"Action 1 (Holonomy) Result: {final_point1}")
    print(f"Action 2 (Deck TF) Result:  {final_point2}")

    if final_point1 == final_point2:
        print("\nThe two actions are the same for the torus T^2.")
    else:
        print("\nThe two actions are different for the torus T^2.")

if __name__ == '__main__':
    solve_actions_on_torus()