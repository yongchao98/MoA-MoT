def solve_torus_action():
    """
    Demonstrates that for the torus T^2, the two actions of the fundamental
    group on the fiber of the universal cover are the same.
    """
    # For the torus X = T^2, its universal cover is the plane, X_tilde = R^2.
    # The fundamental group pi_1(T^2) is isomorphic to Z^2.
    # We represent an element of the fundamental group by a pair of integers (m, n).
    m, n = 2, 3
    gamma = (m, n)

    # The covering map is p(x, y) = (x mod 1, y mod 1).
    # We choose the basepoint x_0 = (0, 0) in T^2.
    # The fiber over x_0 is p^-1(x_0) = Z^2, the set of integer points in R^2.
    # We choose a point x_tilde in the fiber to act upon.
    k, l = 5, 7
    x_tilde = (k, l)

    # We also fix a basepoint in the fiber to define the isomorphism between
    # the fundamental group and the group of deck transformations. Let's use x_tilde_0 = (0, 0).
    x_tilde_0 = (0, 0)

    print(f"Let X = T^2, the 2-torus.")
    print(f"Let the element of the fundamental group be [γ], corresponding to the integer pair (m, n) = {gamma}.")
    print(f"Let the point in the fiber p⁻¹(x₀) be x̃, corresponding to the integer pair (k, l) = {x_tilde}.")
    print("-" * 30)

    # --- Action 1: By restricting deck transformations to the fiber ---
    print("Action 1: By restriction of deck transformations.")
    print("Step 1: Find the deck transformation φ_[γ] corresponding to [γ] = (m, n).")
    # A deck transformation for the cover R^2 -> T^2 is a translation by an integer vector (a, b).
    # The specific transformation φ_[γ] is the one that maps our basepoint x_tilde_0 = (0,0)
    # to the endpoint of the lift of the loop γ starting at x_tilde_0.
    # The lift of γ = (m,n) starting at (0,0) is the path t -> (m*t, n*t) in R^2.
    # Its endpoint at t=1 is (m, n). So, φ_[γ](0,0) = (m, n).
    # A translation by (a, b) maps (0,0) to (a,b).
    # Therefore, the deck transformation is a translation by (m, n).
    deck_translation_vector = gamma
    print(f"The deck transformation φ_[γ] is a translation by the vector (m, n) = {deck_translation_vector}.")

    print("\nStep 2: Apply this deck transformation to the point x̃ = (k, l).")
    result_1_x = x_tilde[0] + deck_translation_vector[0]
    result_1_y = x_tilde[1] + deck_translation_vector[1]
    result_1 = (result_1_x, result_1_y)
    print(f"The resulting point is φ_[γ](x̃) = (k, l) + (m, n) = ({k}, {l}) + ({m}, {n}).")
    print(f"Result 1 = ({k} + {m}, {l} + {n}) = {result_1}.")
    print("-" * 30)

    # --- Action 2: By holonomy around loops in X ---
    print("Action 2: By holonomy (path lifting).")
    print("Step 1: Lift the loop γ, starting from the point x̃ = (k, l).")
    # A loop representing γ is the projection of the path t -> (m*t, n*t).
    # The lift,tilde_γ(t), must satisfy p(tilde_γ(t)) = p(m*t, n*t) and tilde_γ(0) = x̃.
    # The first condition implies tilde_γ(t) = (m*t + c1, n*t + c2) for some integers c1, c2.
    # The second condition (at t=0) implies (c1, c2) = (k, l).
    # So the lifted path is tilde_γ(t) = (m*t + k, n*t + l).
    print(f"The lift of γ starting at ({k}, {l}) is the path t ↦ ({m}*t + {k}, {n}*t + {l}).")

    print("\nStep 2: The action is defined as the endpoint of this lift at t=1.")
    result_2_x = m * 1 + k
    result_2_y = n * 1 + l
    result_2 = (result_2_x, result_2_y)
    print(f"The endpoint at t=1 is ({m}*1 + {k}, {n}*1 + {l}).")
    print(f"Result 2 = ({k} + {m}, {l} + {n}) = {result_2}.")
    print("-" * 30)

    # --- Conclusion ---
    print("Conclusion:")
    if result_1 == result_2:
        print("The results of both actions are identical.")
        print("This holds true for any choice of [γ] in π₁(T²) and any point x̃ in the fiber.")
        print("\nTherefore, the two actions are the same when X = T^2.")
    else:
        # This part should not be reached for the Torus case.
        print("The results are different. The actions are not the same.")

solve_torus_action()