def solve_term_rewriting_completion():
    """
    Solves the Knuth-Bendix completion task.

    This function contains the final list of rules added during the completion
    of the given term-rewriting system. The rules are derived through a process
    of finding critical pairs and inter-reducing the system until it is
    confluent. The final rules are ordered increasingly by their left-hand side
    using the lexicographic path ordering induced by f<g<h.
    """

    # The final set of rules added by the Knuth-Bendix completion algorithm.
    added_rules = [
        "f(f(x,x), f(y,y)) -> f(x,x)",
        "g(x) -> f(x,x)",
        "h(x) -> f(x,x)"
    ]

    # Print the rules, separated by commas, as requested.
    print(", ".join(added_rules))

solve_term_rewriting_completion()