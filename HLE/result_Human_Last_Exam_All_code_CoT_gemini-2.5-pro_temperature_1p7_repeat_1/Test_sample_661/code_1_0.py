def solve_knuth_bendix():
    """
    This function provides the solution to the term-rewriting system completion problem.
    The steps of the Knuth-Bendix algorithm were performed manually as described in the plan.
    This code simply prints the final, derived result.
    """

    # The rules added during the Knuth-Bendix completion process,
    # ordered by their Left-Hand Side (LHS) using the specified LPO (f < g < h).

    # Rule 1 (from reducing original rule R3 with R4)
    rule1_lhs = "f(g(x),g(y))"
    rule1_rhs = "g(x)"

    # Rule 2 (from reducing original rule R1 with R2 and R4)
    rule2_lhs = "g(g(g(x)))"
    rule2_rhs = "g(x)"

    # Rule 3 (from the first critical pair between R1 and R3)
    rule3_lhs = "h(x)"
    rule3_rhs = "g(x)"

    # LPO ordering of LHS: f(g(x),g(y)) < g(g(g(x))) < h(x) due to precedence f < g < h.
    
    final_rules = [
        f"{rule1_lhs} -> {rule1_rhs}",
        f"{rule2_lhs} -> {rule2_rhs}",
        f"{rule3_lhs} -> {rule3_rhs}"
    ]

    print(", ".join(final_rules))

solve_knuth_bendix()
<<<f(g(x),g(y)) -> g(x), g(g(g(x))) -> g(x), h(x) -> g(x)>>>