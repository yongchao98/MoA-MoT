import numpy as np

def solve_torus_actions():
    """
    This function models and compares two actions of the fundamental group
    pi_1(T^2) on the fiber of its universal cover.
    """
    print("This script will demonstrate two actions of the fundamental group on the fiber of a universal cover for X = T^2 (the 2-torus).")
    
    # Mathematical Context
    print("\nLet X be the 2-torus, T^2.")
    print("The universal cover is the plane, tilde{X} = R^2, with the covering map p(x,y) = (x mod 1, y mod 1).")
    print("The fundamental group pi_1(T^2, x_0) is isomorphic to Z^2.")
    print("Let the basepoint x_0 = p(0,0). The fiber p^{-1}(x_0) is the integer lattice Z^2 in R^2.\n")

    # Choose a sample element from the fundamental group and a point in the fiber
    # An element of pi_1(T^2) is represented by an integer vector [m, n].
    loop_class = np.array([3, -2])
    m, n = loop_class

    # A point in the fiber p^{-1}(x_0) is represented by an integer vector [k, l].
    fiber_point = np.array([5, 4])
    k, l = fiber_point

    print(f"We will use a loop class from pi_1(T^2) represented by the vector: [{m}, {n}]")
    print(f"We will see how it acts on a point in the fiber represented by the vector: [{k}, {l}]\n")

    # Action 1: Holonomy around loops (Path Lifting)
    print("--- Action 1: By holonomy (path lifting) ---")
    print("This action is defined by lifting the loop to the universal cover and finding the endpoint.")
    print(f"A representative loop in T^2 is gamma(t) = p({m}*t, {n}*t) for t in [0,1].")
    print(f"We lift this path to R^2 starting at the fiber point ({k}, {l}).")
    # The lifted path is tilde_gamma(t) = (m*t + k, n*t + l).
    print(f"The unique lifted path tilde_gamma(t) is: ({m}*t + {k}, {n}*t + {l}).")
    # The action is defined as the endpoint of this path at t=1.
    result_action1 = fiber_point + loop_class
    print(f"The endpoint at t=1 is tilde_gamma(1) = ({m}*1 + {k}, {n}*1 + {l}) = ({result_action1[0]}, {result_action1[1]}).")
    print(f"Result of Action 1: [{result_action1[0]}, {result_action1[1]}]\n")


    # Action 2: Restricting deck transformations to the fiber
    print("--- Action 2: By restricting deck transformations ---")
    print("This action is defined by finding the deck transformation for the loop class, then applying it.")

    # A deck transformation for the torus is a translation by an integer vector.
    # To find the right transformation for the loop [m, n], we lift the loop from a base fiber point, typically (0,0).
    base_fiber_point = np.array([0, 0])
    print(f"To identify the deck transformation phi, we see where it maps the base fiber point ({base_fiber_point[0]},{base_fiber_point[1]}).")
    print(f"This is determined by the endpoint of the loop [{m},{n}] when lifted starting at ({base_fiber_point[0]},{base_fiber_point[1]}).")
    
    # The endpoint of the lift from (0,0) is simply (m,n).
    translation_vector = loop_class
    print(f"The lift ends at ({translation_vector[0]}, {translation_vector[1]}). Thus, the corresponding deck transformation phi is a translation by this vector.")
    print(f"phi(x, y) = (x + {m}, y + {n}).")

    # Now, we apply this deck transformation to our fiber point (k, l).
    result_action2 = fiber_point + translation_vector
    print(f"Applying phi to our fiber point ({k}, {l}):")
    print(f"phi({k}, {l}) = ({k} + {m}, {l} + {n}) = ({result_action2[0]}, {result_action2[1]}).")
    print(f"Result of Action 2: [{result_action2[0]}, {result_action2[1]}]\n")

    # Compare the results
    print("--- Conclusion ---")
    print(f"Result from Action 1 (Holonomy): [{result_action1[0]}, {result_action1[1]}]")
    print(f"Result from Action 2 (Deck Transformations): [{result_action2[0]}, {result_action2[1]}]")

    if np.array_equal(result_action1, result_action2):
        print("\nThe two actions produce the exact same result.")
    else:
        print("\nThe two actions produce different results.")

solve_torus_actions()