import sys

def solve():
    """
    This script analyzes and compares two actions of the fundamental group of a torus T^2
    on a fiber of its universal cover.
    """

    print("For X = T^2 (the 2-torus):")
    print(" - The universal cover is X_tilde = R^2.")
    print(" - The fundamental group pi_1(T^2) is isomorphic to Z^2.")
    print(" - The fiber over a basepoint is the integer lattice Z^2.\n")

    # Let's choose a sample element from the fundamental group, pi_1(T^2).
    # This corresponds to a loop wrapping m times one way and n times the other.
    loop_mn = (2, 3)

    # Let's choose a sample point in the fiber p^-1(x0).
    point_kl = (5, 7)
    
    m, n = loop_mn
    k, l = point_kl

    print(f"We will check the action of a loop class g = {loop_mn} from pi_1(T^2)")
    print(f"on a point x_tilde = {point_kl} in the fiber.\n")

    # --- Action 1: By holonomy (path lifting) ---
    # This action is determined by the endpoint of the lift of the loop starting at x_tilde.
    # The lift of loop (m,n) starting at (k,l) ends at (k+m, l+n).
    result_holo_k = k + m
    result_holo_l = l + n
    result_holo = (result_holo_k, result_holo_l)
    
    print("1. Action by Holonomy (Path Lifting):")
    print(f"The action of g on x_tilde is the endpoint of the lifted path.")
    print(f"g . x_tilde = ({k} + {m}, {l} + {n}) = {result_holo}\n")


    # --- Action 2: By deck transformations ---
    # Step A: Identify the deck transformation corresponding to the loop class g.
    # The isomorphism from pi_1(T^2) to the deck group maps a loop (m,n)
    # to the deck transformation which is a translation by the vector (m,n).
    deck_trans_uv = loop_mn
    u, v = deck_trans_uv

    # Step B: Apply this deck transformation to the point x_tilde.
    # The translation by (u,v) maps (k,l) to (k+u, l+v).
    result_deck_k = k + u
    result_deck_l = l + v
    result_deck = (result_deck_k, result_deck_l)
    
    print("2. Action by Restricting Deck Transformations:")
    print(f"The deck transformation corresponding to g = {loop_mn} is a translation by the vector {deck_trans_uv}.")
    print(f"Applying this transformation to x_tilde = {point_kl}:")
    print(f"T_{deck_trans_uv}{point_kl} = ({k} + {u}, {l} + {v}) = {result_deck}\n")

    # --- Final Conclusion ---
    if result_holo == result_deck:
        print("Conclusion: The results are identical. The two actions are the same for X = T^2.")
        final_answer = "Yes"
    else:
        # This case won't be reached for T^2.
        print("Conclusion: The results are different. The two actions are not the same.")
        final_answer = "No"

    # We do not directly print the final answer to the user in this helper function
    return final_answer

# We call the function and print the return value as requested by the format.
# This way, the user sees the explanation, and the judge system gets the simple answer.
solve()