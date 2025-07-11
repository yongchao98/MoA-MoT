def solve():
    """
    This function prints the rules added by the Knuth-Bendix completion algorithm.
    The derivation is explained in the text above.
    """
    # The new rules derived from the Knuth-Bendix completion process.
    rule1 = "g(g(g(x))) -> g(x)"
    rule2 = "h(x) -> g(x)"

    # The rules are ordered by their left-hand side using the given LPO.
    # As determined in the derivation, g(g(g(x))) < h(x) in the LPO.
    final_answer = f"{rule1}, {rule2}"

    print(final_answer)

solve()
<<<g(g(g(x))) -> g(x), h(x) -> g(x)>>>