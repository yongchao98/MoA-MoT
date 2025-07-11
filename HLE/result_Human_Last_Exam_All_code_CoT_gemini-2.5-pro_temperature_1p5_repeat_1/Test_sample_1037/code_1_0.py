def solve_scale_cardinalities():
    """
    This function calculates and prints the cardinalities based on the analysis of the category of scales.
    """
    
    # 1. Cardinality of S/A
    # The quotient space S/A is the group *R/Z.
    # The cardinality of the hyperreals, |*R|, is 2^c, where c = 2^aleph_0 = Beth_1.
    # So, |*R| = 2^(Beth_1) = Beth_2.
    # Since |*R| = |*R/Z| * |Z|, we have Beth_2 = |*R/Z| * Beth_0, which implies |*R/Z| = Beth_2.
    card_S_A = "Beth_2"
    
    # 2. Cardinality of B/S
    # The terminal object B is the trivial scale Z -> {0}. Its group is {0}.
    # The quotient B/S is {0} / im(h: *R -> {0}), which is {0}/{0} = {0}.
    # The cardinality of the trivial group {0} is 1.
    card_B_S = 1
    
    # 3. Cardinality of H_1(B/A, Q)
    # The quotient space B/A is {0} / im(h: Z -> {0}), which is {0}/{0}.
    # This space is a single point.
    # The first homology group of a point, H_1({0}, Q), is the trivial group 0.
    # The cardinality of this group is 1.
    card_H1_B_A = 1
    
    # The final answer is the sequence of these three cardinalities.
    print(f"{card_S_A} {card_B_S} {card_H1_B_A}")

solve_scale_cardinalities()