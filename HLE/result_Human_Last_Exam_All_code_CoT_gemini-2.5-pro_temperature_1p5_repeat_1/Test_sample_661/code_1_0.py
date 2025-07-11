def solve():
    """
    This function prints the solution to the term-rewriting system problem.
    The rules are derived from the Knuth-Bendix completion algorithm and ordered
    by their left-hand sides using the specified lexicographic path ordering.
    """

    # The rules added by the Knuth-Bendix completion, ordered as determined.
    added_rules = [
        "f(y, y) -> g(g(y))",
        "f(g(x), g(y)) -> g(x)",
        "g(g(g(x))) -> g(x)",
        "h(x) -> g(x)"
    ]

    # Joining the rules with a comma and a space for the final output.
    result = ", ".join(added_rules)
    print(result)

solve()