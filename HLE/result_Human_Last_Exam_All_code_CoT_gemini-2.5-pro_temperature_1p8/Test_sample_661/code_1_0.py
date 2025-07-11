def solve():
    """
    This function prints the new rules added by the Knuth-Bendix completion process,
    ordered increasingly by their left-hand side according to the specified term ordering.
    """
    # The derived new rules are:
    # 1. f(g(x), g(y)) -> g(x)
    # 2. g(g(y)) -> f(y, y)
    # 3. h(x) -> g(x)
    
    # We format them as requested.
    rule1 = "f(g(x), g(y)) -> g(x)"
    rule2 = "g(g(y)) -> f(y, y)"
    rule3 = "h(x) -> g(x)"
    
    result = [rule1, rule2, rule3]
    
    print(", ".join(result))

solve()
<<<f(g(x), g(y)) -> g(x), g(g(y)) -> f(y, y), h(x) -> g(x)>>>