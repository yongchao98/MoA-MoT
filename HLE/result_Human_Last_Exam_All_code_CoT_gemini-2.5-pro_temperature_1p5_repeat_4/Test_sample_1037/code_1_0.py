def solve_cardinalities():
    """
    Solves for the cardinalities of S/A, B/S, and H_1(B/A, Q).

    The solution proceeds as follows:
    1.  The initial object A is the scale (id: Z -> Z). Its group G_A has cardinality Beth_0.
    2.  The scale S is the inclusion (incl: Z -> *R). Its group G_S has cardinality Beth_1.
    3.  The quotient S/A is G_S / im(h_AS) = *R / Z. Its cardinality is |*R| = Beth_1.
    4.  The terminal object B must be universal. Following the pattern of cardinalities from A (Beth_0)
        and S (Beth_1), its group G_B is inferred to have cardinality Beth_2.
    5.  The quotient B/S is G_B / im(h_SB). Since |im(h_SB)| <= |*R| = Beth_1 and |G_B| = Beth_2,
        the cardinality of the quotient |B/S| is Beth_2.
    6.  The quantity H_1(B/A, Q) is the torsion-free rank of G_B / im(h_AB) = G_B / Z.
        This rank is Rank(G_B) - Rank(Z) = Beth_2 - 1 = Beth_2.
    """

    # Value for |S/A|
    card_S_div_A = "Beth_1"

    # Value for |B/S|
    card_B_div_S = "Beth_2"

    # Value for H_1(B/A, Q)
    card_H1_B_div_A = "Beth_2"

    print(f"{card_S_div_A} {card_B_div_S} {card_H1_B_div_A}")

solve_cardinalities()