import numpy as np

def solve_torus_action():
    """
    Demonstrates that for X = T^2, the action of pi_1 on a fiber
    by holonomy is the same as the action by deck transformations.
    """
    # Let X = T^2 be the 2-torus.
    # The universal cover X_tilde is R^2.
    # The base point x_0 on T^2 is the image of the origin in R^2.
    # The fiber p^-1(x_0) is the integer lattice Z^2.
    # The fundamental group pi_1(T^2, x_0) is isomorphic to Z^2.

    # We choose a sample element of the fundamental group, represented by an integer vector.
    # This corresponds to a loop wrapping k times one way and l times the other.
    gamma_class = np.array([2, 3])

    # We choose a sample point in the fiber p^-1(x_0).
    fiber_point = np.array([5, 7])
    
    # We choose a base point in the fiber to define the isomorphism between pi_1 and the deck group.
    fiber_base_point = np.array([0, 0])

    print(f"We test the two actions for X = T^2.")
    print(f"Chosen fundamental group element [gamma] ~ {list(gamma_class)}")
    print(f"Chosen fiber point ~ {list(fiber_point)}\n")

    # --- Action 1: Holonomy (Path Lifting) ---

    # For T^2, lifting the loop corresponding to vector [k, l] starting from [m, n]
    # results in a path from [m, n] to [m+k, n+l].
    holonomy_result = fiber_point + gamma_class
    
    k, l = gamma_class
    m, n = fiber_point
    
    print("Action 1 (Holonomy):")
    print("The action is defined by the endpoint of the lifted path.")
    print(f"Lifting the loop [gamma] starting at ({m}, {n}) results in a path ending at:")
    print(f"({m} + {k}, {n} + {l}) = {list(holonomy_result)}\n")


    # --- Action 2: Deck Transformations ---

    # 1. Find the deck transformation g_gamma corresponding to [gamma].
    #    g_gamma is the unique deck transformation mapping fiber_base_point to the
    #    endpoint of the lift of gamma starting at fiber_base_point.
    #    The deck group for T^2 consists of translations by integer vectors.
    lift_endpoint_from_base = fiber_base_point + gamma_class
    
    # The deck transformation g_gamma is thus translation by the vector `lift_endpoint_from_base`.
    # Since fiber_base_point is (0,0), this is just translation by gamma_class.
    translation_vector = gamma_class
    
    # 2. Apply this deck transformation to the fiber_point.
    deck_action_result = fiber_point + translation_vector
    
    print("Action 2 (Deck Transformations):")
    print("The action is defined by first finding the corresponding deck transformation, then applying it.")
    print(f"1. The deck transformation g_gamma for [gamma] ~ {list(gamma_class)} is translation by the vector {list(translation_vector)}.")
    print(f"2. Applying this translation to the fiber point ({m}, {n}):")
    print(f"({m}, {n}) + ({k}, {l}) = {list(deck_action_result)}\n")


    # --- Comparison ---
    print("--- Comparison ---")
    print(f"Result of Action 1 (Holonomy): {list(holonomy_result)}")
    print(f"Result of Action 2 (Deck Transformations): {list(deck_action_result)}")

    are_same = np.array_equal(holonomy_result, deck_action_result)

    if are_same:
        print("\nThe two actions are the same.")
    else:
        print("\nThe two actions are NOT the same.")

solve_torus_action()