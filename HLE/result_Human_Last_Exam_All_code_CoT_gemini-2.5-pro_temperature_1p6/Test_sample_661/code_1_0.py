def solve_term_rewriting():
    """
    This function prints the ordered list of rules added by the Knuth-Bendix completion.
    The derivation is based on the logic described in the plan.
    """
    # The new rules derived from the Knuth-Bendix completion procedure.
    rule1_lhs = "f(g(x), g(y))"
    rule1_rhs = "g(x)"

    rule2_lhs = "g(g(g(x)))"
    rule2_rhs = "g(x)"

    rule3_lhs = "h(x)"
    rule3_rhs = "g(x)"
    
    # Based on the LPO f < g < h, the LHS are ordered as:
    # f(g(x), g(y)) < g(g(g(x))) < h(x)
    ordered_rules = [
        f"{rule1_lhs} -> {rule1_rhs}",
        f"{rule2_lhs} -> {rule2_rhs}",
        f"{rule3_lhs} -> {rule3_rhs}"
    ]

    print(", ".join(ordered_rules))

solve_term_rewriting()