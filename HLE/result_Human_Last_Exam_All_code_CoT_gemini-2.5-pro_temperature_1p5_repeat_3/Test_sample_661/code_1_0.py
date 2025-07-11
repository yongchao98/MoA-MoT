def knuth_bendix_completion_result():
    """
    This function presents the results of the Knuth-Bendix completion algorithm
    for the given term-rewriting system and ordering.
    """

    # The completion process generates new rules from critical pairs and by simplifying
    # existing rules. The following rules are the result of this process.

    # 1. From the critical pair between f(g(x),h(x))->g(x) and f(g(x),h(y))->h(x),
    #    we deduce g(x)=h(x). With the ordering f<g<h, this gives:
    rule_A = "h(x) -> g(x)"

    # 2. Simplifying the original rules with rule_A yields the next two rules.
    #    f(y,y)->g(h(y)) becomes:
    rule_B = "f(y, y) -> g(g(y))"
    #    f(g(x),h(y))->h(x) becomes:
    rule_C = "f(g(x), g(y)) -> g(x)"

    # 3. A new critical pair between rule_B and rule_C (on the term f(g(z),g(z)))
    #    yields g(g(g(z)))=g(z). This orients to:
    rule_D = "g(g(g(x))) -> g(x)"

    # The problem asks to list all added rules, ordered by their left-hand side (LHS)
    # using the lexicographic path ordering (LPO) with f < g < h.
    # The LHSs are:
    # - f(y, y)
    # - f(g(x), g(y))
    # - g(g(g(x)))
    # - h(x)
    #
    # The LPO gives the following ascending order for the LHSs:
    # f(y, y) < f(g(x), g(y)) < g(g(g(x))) < h(x)
    #
    # So, we list the rules in this order.
    
    final_ordered_rules = [
        rule_B,
        rule_C,
        rule_D,
        rule_A
    ]

    # Print the rules, separated by commas.
    print(", ".join(final_ordered_rules))

knuth_bendix_completion_result()