import sys

def solve():
    """
    This script determines if the holonomy action and the deck transformation action
    of the fundamental group on the fiber of the universal cover are the same for X = T^2.
    It uses a concrete example to demonstrate the calculation for both actions.
    """

    # --- Setup ---
    # We choose a specific element of the fundamental group and a specific point in the fiber
    # to make the demonstration concrete.

    # An element of pi_1(T^2, x_0) can be represented by an integer pair (m, n).
    m, n = 3, 5
    gamma_loop_class = (m, n)

    # A point in the fiber p^-1(x_0) is a point in the integer lattice Z^2.
    k, l = 2, 7
    fiber_point = (k, l)
    
    # Base point in the universal cover R^2
    x0_tilde = (0, 0)

    print("We investigate the two actions for the case X = T^2, the 2-torus.")
    print("The universal cover is X_tilde = R^2.")
    print("The fundamental group pi_1(T^2) is isomorphic to Z^2.")
    print("The fiber p^-1(x_0) over the basepoint is the integer lattice Z^2 in R^2.")
    print(f"Let's choose a loop class [gamma] corresponding to the integer vector ({m}, {n}).")
    print(f"Let's choose a fiber point y_tilde = ({k}, {l}).\n")

    # --- Action 1: Holonomy ---
    print("--- 1. Action by Holonomy ---")
    print("The action of [gamma] on y_tilde is the endpoint of the lift of the loop gamma that starts at y_tilde.")
    print(f"A representative loop gamma(t) in T^2 can be visualized as the projection of a straight path from (0,0) to ({m},{n}) in R^2.")
    print(f"The unique lift of this loop to R^2 that starts at y_tilde = ({k},{l}) is the pathtilde_gamma(t) = ({k} + {m}*t, {l} + {n}*t).")
    print(f"To find the result of the action, we evaluate the endpoint of this path at t=1:")
    print(f"Endpoint = ({k} + {m}, {l} + {n})")
    
    holonomy_result_k = k + m
    holonomy_result_l = l + n
    print(f"Resulting Point = ({holonomy_result_k}, {holonomy_result_l})\n")

    # --- Action 2: Deck Transformations ---
    print("--- 2. Action by Restricting Deck Transformations ---")
    print("First, we determine the deck transformation phi_[gamma] that corresponds to the loop class [gamma].")
    print(f"By definition, phi_[gamma] is the unique deck transformation that maps the basepoint of the cover, x0_tilde = {x0_tilde}, to the endpoint of the lift of gamma starting at x0_tilde.")
    print(f"The lift of gamma starting at (0,0) is the path ({m}*t, {n}*t). Its endpoint at t=1 is ({m}, {n}).")
    print(f"Deck transformations for the torus are translations by integer vectors, T_(a,b)(x,y) = (x+a, y+b).")
    print(f"The transformation that maps (0,0) to ({m},{n}) is T_({m},{n}). So, phi_[gamma] = T_({m},{n}).")
    print(f"Now, we apply this deck transformation to our fiber point y_tilde = ({k},{l}):")
    print(f"Action = T_({m},{n})(({k},{l})) = ({k} + {m}, {l} + {n})")

    deck_result_k = k + m
    deck_result_l = l + n
    print(f"Resulting Point = ({deck_result_k}, {deck_result_l})\n")

    # --- Conclusion ---
    print("--- Conclusion ---")
    print(f"Holonomy action result: ({holonomy_result_k}, {holonomy_result_l})")
    print(f"Deck transformation action result: ({deck_result_k}, {deck_result_l})")

    if (holonomy_result_k, holonomy_result_l) == (deck_result_k, deck_result_l):
        print("The results are identical. The two actions are the same for X = T^2.")
    else:
        print("The results are different. The two actions are not the same for X = T^2.")

solve()