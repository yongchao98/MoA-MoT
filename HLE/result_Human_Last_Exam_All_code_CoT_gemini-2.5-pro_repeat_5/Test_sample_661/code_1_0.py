def solve_knuth_bendix():
    """
    This function prints the rules added by the Knuth-Bendix completion algorithm.
    The derivation is explained in the text preceding this code block.
    """
    
    # The final set of rules derived from the completion process
    rules = [
        "f(y, y) -> g(g(y))",
        "f(g(x), g(y)) -> g(x)",
        "g(g(g(x))) -> g(x)",
        "h(x) -> g(x)"
    ]
    
    # The rules are already ordered by their left-hand side using the LPO f<g<h.
    # We join them with commas as requested.
    result = ", ".join(rules)
    
    print(result)

solve_knuth_bendix()
<<<f(y, y) -> g(g(y)), f(g(x), g(y)) -> g(x), g(g(g(x))) -> g(x), h(x) -> g(x)>>>