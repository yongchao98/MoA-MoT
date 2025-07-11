def solve_category_theory_problem():
    """
    This function calculates the cardinalities of three group quotients
    derived from objects in the category of scales.

    The final equation of cardinalities is |S/A|, |B/S|, and |H_1(B/A, Q)|.
    """

    # Part 1: Cardinality of S/A
    # A is the initial object, the scale id: Z -> Z. So its space is Z.
    # S is the scale of inclusion Z -> *R (the hyperreals). Its space is *R.
    # The canonical map A -> S is the inclusion Z -> *R.
    # The quotient S/A is the space *R / Z.
    # The cardinality of Z is aleph_0, or Beth_0.
    # The cardinality of the hyperreals *R (via standard ultraproduct construction)
    # is the cardinality of the continuum, c = 2^aleph_0, which is Beth_1.
    # Using cardinal arithmetic, |*R| = |*R/Z| * |Z|.
    # This means Beth_1 = max(|*R/Z|, Beth_0), which implies |*R/Z| = Beth_1.
    card_S_div_A = "Beth_1"

    # Part 2: Cardinality of B/S
    # B is the terminal object, the scale zero: Z -> {0}. Its space is {0}.
    # S has the space *R.
    # The canonical map S -> B is the zero map *R -> {0}. Its image is {0}.
    # The quotient B/S is the space {0} / {0}.
    # This is the trivial group, which contains one element.
    card_B_div_S = 1

    # Part 3: Cardinality of H_1(B/A, Q)
    # B's space is {0} and A's space is Z.
    # The canonical map A -> B is the zero map Z -> {0}. Its image is {0}.
    # The quotient space B/A is {0} / {0}, which is a single point space.
    # We need the cardinality of the first singular homology group of a point
    # with coefficients in the rational numbers Q, H_1({pt}, Q).
    # The homology of a point is trivial for dimensions n > 0.
    # So, H_1({pt}, Q) is the trivial group {0}.
    # The cardinality of the trivial group is 1.
    card_H1_B_div_A = 1

    # The final answer expressed in the required format.
    print(f"{card_S_div_A} {card_B_div_S} {card_H1_B_div_A}")

solve_category_theory_problem()