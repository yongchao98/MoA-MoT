def solve_knuth_bendix():
    """
    This function provides the solution to the Knuth-Bendix completion problem.
    It constructs the final list of added rules and prints them in the specified format.
    """
    # The new rules generated and simplified by the Knuth-Bendix completion algorithm are:
    # 1. f(y, y) -> g(g(y))
    # 2. f(g(x), g(y)) -> g(x)
    # 3. g(g(g(x))) -> g(x)
    # 4. h(x) -> g(x)
    # These have been ordered by their left-hand side according to the LPO f<g<h.

    rules = [
        "f(y, y) -> g(g(y))",
        "f(g(x), g(y)) -> g(x)",
        "g(g(g(x))) -> g(x)",
        "h(x) -> g(x)"
    ]

    # The final answer is the comma-separated list of these rules.
    final_answer = ", ".join(rules)
    print(final_answer)

solve_knuth_bendix()