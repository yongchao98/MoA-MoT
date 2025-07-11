import sympy

def explain_actions_on_torus():
    """
    Analyzes and compares two actions of the fundamental group of a torus
    on the fiber of its universal cover.
    """

    # Represent elements symbolically
    m, n, k, l = sympy.symbols('m n k l', integer=True)
    
    # Let X be the 2-torus T^2. Its universal cover is R^2.
    # The fundamental group pi_1(T^2) is isomorphic to Z^2.
    # An element of pi_1(T^2) is represented by a vector (m, n).
    gamma = sympy.Matrix([m, n])

    # The fiber p^{-1}(x_0) over the basepoint is the integer lattice Z^2.
    # A point in the fiber is represented by a vector (k, l).
    x_tilde = sympy.Matrix([k, l])

    print("Let X be the 2-torus T^2.")
    print("The fundamental group pi_1(T^2) is isomorphic to the group of integer vectors Z^2.")
    print(f"Let an element of the fundamental group be represented by the vector g = {gamma.T}.")
    print("\nThe fiber of the universal cover over the basepoint is the integer lattice Z^2.")
    print(f"Let a point in the fiber be represented by the vector x = {x_tilde.T}.")
    print("-" * 30)

    # --- Action 1: Holonomy ---
    print("\nAction 1: Action by holonomy around loops.")
    print("This action is defined by lifting the loop corresponding to g starting from the point x.")
    print("A path in the universal cover R^2 corresponding to the loop g is p(t) = t*g.")
    print("The lift of this loop starting at x is p_lift(t) = x + t*g.")
    print(f"p_lift(t) = {x_tilde.T} + t*{gamma.T} = {x_tilde.T + sympy.Symbol('t')*gamma.T}")
    print("The endpoint of this lifted path is at t=1.")
    holonomy_result = x_tilde + 1 * gamma
    print(f"Endpoint = x + g = {x_tilde.T} + {gamma.T} = {holonomy_result.T}")
    
    # --- Action 2: Deck Transformations ---
    print("\nAction 2: Action by restricting deck transformations to the fiber.")
    print("For the torus T^2, the deck transformations are translations by integer vectors.")
    print("There is a canonical isomorphism between pi_1(T^2) and the group of deck transformations.")
    print("The group element g in pi_1(T^2) corresponds to the deck transformation T_g, which is a translation by the vector g.")
    print(f"The action of g on x is defined by applying the deck transformation T_g to x.")
    deck_result = x_tilde + gamma # T_g(x) = x + g
    print(f"T_g(x) = x + g = {x_tilde.T} + {gamma.T} = {deck_result.T}")
    
    print("-" * 30)
    print("\nComparing the results of the two actions:")
    print(f"Action 1 (Holonomy) result: {holonomy_result.T}")
    print(f"Action 2 (Deck) result:    {deck_result.T}")
    
    are_same = (holonomy_result == deck_result)
    
    print(f"\nThe two actions produce the same result.")
    print("Therefore, for X = T^2, the two actions are the same.")
    
    if are_same:
        answer = "Yes"
    else:
        answer = "No"
    
    print(f"\nFinal Answer: {answer}")


explain_actions_on_torus()
<<<Yes>>>