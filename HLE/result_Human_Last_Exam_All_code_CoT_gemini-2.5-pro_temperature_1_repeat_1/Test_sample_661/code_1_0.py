def solve_knuth_bendix_completion():
    """
    This function prints the result of the Knuth-Bendix completion.
    The completion process was performed manually based on the given term-rewriting
    system and term ordering. The added rules are ordered according to the
    lexicographic path ordering (LPO) with precedence f < g < h.
    """

    # The new rules derived from the completion algorithm.
    # The order is determined by comparing the left-hand sides (LHS)
    # using the specified LPO. The ordering of the LHSs is:
    # f(f(x,x), f(x,x)) < f(g(x), g(y)) < g(f(x,x)) < g(g(y)) < h(x)
    added_rules = [
        "f(f(x,x), f(x,x)) -> f(x,x)",
        "f(g(x), g(y)) -> g(x)",
        "g(f(x,x)) -> g(x)",
        "g(g(y)) -> f(y,y)",
        "h(x) -> g(x)"
    ]

    # Joining the rules into a single string, separated by commas.
    result = ", ".join(added_rules)
    print(result)

solve_knuth_bendix_completion()