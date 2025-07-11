def solve():
    """
    Calculates the satisfaction scores s(N,W1) and s(N,W2) and their ratio.
    """

    # Based on the step-by-step derivation:
    # W1 is a committee in the core with the lowest satisfaction for group N.
    # The core property imposes very strong stability conditions. To prevent
    # various subgroups of N from defecting, each member i of N must have
    # a satisfaction of at least 8.
    # The minimum satisfaction for group N under this constraint is achieved
    # with a committee where 7 common candidates are chosen, and 1 specific
    # candidate for each of the 8 members of N is chosen.
    # s(N,W1) = 8 * w_c + sum(w_i) = 8 * 7 + 8 * 1 = 56 + 8 = 64.
    s_N_W1 = 64

    # W2 is a committee satisfying Extended Justified Representation (EJR)
    # with the lowest satisfaction for group N.
    # EJR for the 8-cohesive group N and its subgroups of size >= 4 implies
    # that at least 5 members of N must have a satisfaction of at least 8.
    # The satisfaction for N is minimized by picking 6 common candidates, and
    # giving 2 specific candidates to each of these 5 members of N.
    # The remaining 3 members of N get a satisfaction of 6.
    # s(N,W2) = 8 * w_c + sum(w_i) = 8 * 6 + 5 * 2 = 48 + 10 = 58.
    s_N_W2 = 58

    # The ratio is s(N,W1) / s(N,W2).
    ratio = s_N_W1 / s_N_W2

    print(f"The lowest satisfaction for N with a core committee is s(N,W1) = {s_N_W1}.")
    print(f"The lowest satisfaction for N with an EJR committee is s(N,W2) = {s_N_W2}.")
    print(f"The ratio is s(N,W1) / s(N,W2) = {s_N_W1} / {s_N_W2} = {ratio}")
    # The final answer in the required format.
    print(f"\n<<<{ratio}>>>")

solve()