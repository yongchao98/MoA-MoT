def knuth_bendix_completion():
    """
    This function follows the steps of the Knuth-Bendix completion algorithm
    for the specific problem provided and prints the final list of added rules.
    """

    # The problem defines three initial rules and a term ordering.
    # R1: f(g(x), h(x)) -> g(x)
    # R2: f(y, y) -> g(h(y))
    # R3: f(g(x), h(y)) -> h(x)
    # Term ordering is LPO induced by f < g < h.

    # Step 1: Find the first critical pair from the initial rules.
    # A critical pair arises from overlapping the left-hand sides (LHS) of R1 and R3.
    # Unifying f(g(x), h(x)) and f(g(x'), h(y')) leads to the critical pair (g(x), h(x)).
    # Using the LPO with h > g, we orient this pair into a new rule.
    new_rule_1 = "h(x) -> g(x)"

    # Step 2: Simplify (inter-reduce) the existing rules with the new rule.
    # The new rule h(x) -> g(x) simplifies R2 and R3.

    # R2, f(y, y) -> g(h(y)), is reduced on its right-hand side.
    # g(h(y)) becomes g(g(y)). The modified rule is:
    rule_from_R2_reduction = "f(y, y) -> g(g(y))"

    # R3, f(g(x), h(y)) -> h(x), is reduced on both sides.
    # The LHS f(g(x), h(y)) becomes f(g(x), g(y)).
    # The RHS h(x) becomes g(x). The modified rule is:
    rule_from_R3_reduction = "f(g(x), g(y)) -> g(x)"

    # The original R1, f(g(x), h(x)) -> g(x), is also reduced to f(g(x), g(x)) -> g(x),
    # but this is a special case of the reduced R3, so it is redundant.

    # Step 3: Find the next critical pair from the current rule set.
    # A new critical pair is found by overlapping the two rules derived above:
    # - f(y, y) -> g(g(y))
    # - f(g(x), g(y)) -> g(x)
    # Overlapping them on the term f(g(x), g(x)) gives the pair (g(g(g(x))), g(x)).
    # Using LPO, g(g(g(x))) > g(x). This orients into our final new rule.
    new_rule_2 = "g(g(g(x))) -> g(x)"

    # Step 4: Finalize and order the added rules.
    # The process is now complete as no new critical pairs can be generated.
    # The set of all added rules is:
    # - h(x) -> g(x)
    # - f(y, y) -> g(g(y))
    # - f(g(x), g(y)) -> g(x)
    # - g(g(g(x))) -> g(x)
    # We order them by their LHS using the LPO f < g < h.
    # The sorted order of LHSs is: f(y, y) < f(g(x), g(y)) < g(g(g(x))) < h(x).
    
    sorted_added_rules = [
        rule_from_R2_reduction,
        rule_from_R3_reduction,
        new_rule_2,
        new_rule_1
    ]

    print(", ".join(sorted_added_rules))

knuth_bendix_completion()