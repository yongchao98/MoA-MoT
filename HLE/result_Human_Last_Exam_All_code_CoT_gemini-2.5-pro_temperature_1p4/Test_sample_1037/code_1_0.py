def solve_scale_cardinalities():
    """
    This function determines the cardinalities based on the mathematical reasoning outlined above.
    The calculations involve concepts from category theory, group theory, and algebraic topology.
    """

    # Part 1: Cardinality of S/A
    # S/A is the quotient group *R / Z.
    # |*R| = Beth_1 and |Z| = Beth_0.
    # From |*R/Z| * |Z| = |*R|, it follows that |S/A| = Beth_1.
    card_S_div_A = "Beth_1"

    # Part 2: Cardinality of B/S
    # B/S is interpreted as the mapping cone of the map f: *R -> {0}.
    # The cardinality of this space is |*R x (0,1]| + 1, which is Beth_1.
    card_B_div_S = "Beth_1"

    # Part 3: Cardinality of H_1(B/A, Q)
    # B/A is the mapping cone of the map f: Z -> {0}.
    # H_1(B/A, Q) is a Q-vector space of countable dimension.
    # Such a vector space has cardinality aleph_0, which is Beth_0.
    card_H1_B_div_A = "Beth_0"

    # Print the final result in the required format
    print(f"{card_S_div_A} {card_B_div_S} {card_H1_B_div_A}")

solve_scale_cardinalities()