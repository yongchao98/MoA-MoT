def print_knuth_bendix_rules():
    """
    This function stores and prints the final set of rules derived from the
    Knuth-Bendix completion algorithm, ordered as requested.
    """
    # The set of new rules derived from the completion process.
    # These are pre-calculated and ordered according to the LPO f < g < h.
    rules = [
        "f(f(x, x), f(y, y)) -> f(x, x)",
        "f(f(y, y), g(z)) -> f(y, y)",
        "f(g(x), f(y, y)) -> g(x)",
        "f(g(x), g(y)) -> g(x)",
        "g(g(y)) -> f(y, y)",
        "h(x) -> g(x)"
    ]

    # Print the rules separated by a comma.
    print(", ".join(rules))

print_knuth_bendix_rules()
<<<f(f(x, x), f(y, y)) -> f(x, x), f(f(y, y), g(z)) -> f(y, y), f(g(x), f(y, y)) -> g(x), f(g(x), g(y)) -> g(x), g(g(y)) -> f(y, y), h(x) -> g(x)>>>