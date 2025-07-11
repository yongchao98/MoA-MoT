def solve_math_problem():
    """
    This function calculates the cardinalities requested in the problem.

    Step-by-step derivation:
    1.  Initial Object A: The initial scale is A = (Z, id_Z), where id_Z is the identity map on the integers.
    2.  Terminal Object B: The category of scales (nontrivial maps) does not have a terminal object. We assume the problem implies a terminal object B in the larger category of all maps, which is the trivial map f_B: Z -> {0}.
    3.  Scale S: S is the inclusion map from Z into the hyperreals, *R. The codomain is *R.
    4.  Canonical Maps and Quotients:
        - S/A: The map A -> S is the inclusion Z -> *R. The quotient is *R / Z.
        - B/S: The map S -> B is the zero map *R -> {0}. The quotient is {0} / {0}.
        - B/A: The map A -> B is the zero map Z -> {0}. The quotient is {0} / {0}.
    5.  Cardinality Calculations:
        - |S/A| = |*R / Z|: The cardinality of the hyperreals *R is the same as the continuum, c = 2^aleph_0 = Beth_1. The cardinality of the quotient *R / Z is also Beth_1.
        - |B/S| = |{0} / {0}|: This is the trivial group, with cardinality 1.
        - |H_1(B/A, Q)|: The space B/A is a single point. The first homology group of a point, H_1(pt, Q), is the trivial group {0}, which has cardinality 1.
    """

    # Cardinality of S/A is Beth_1
    card_S_A = "Beth_1"

    # Cardinality of B/S is 1
    card_B_S = 1

    # Cardinality of H_1(B/A, Q) is 1
    card_H1_B_A = 1

    # Print the final answer in the required format
    print(f"{card_S_A} {card_B_S} {card_H1_B_A}")

solve_math_problem()