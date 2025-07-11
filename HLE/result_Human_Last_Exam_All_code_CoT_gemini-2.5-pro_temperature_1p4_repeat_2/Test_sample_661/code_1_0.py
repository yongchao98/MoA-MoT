def solve_knuth_bendix():
    """
    This function solves the Knuth-Bendix completion problem as described.
    The process involves finding critical pairs, simplifying and orienting them into new rules,
    and then ordering the final set of added rules.
    The derivation steps are as follows:

    1. Initial Rules:
       R1: f(g(x), h(x)) -> g(x)
       R2: f(y, y) -> g(h(y))
       R3: f(g(x), h(y)) -> h(x)
       Ordering: LPO with f < g < h

    2. First Critical Pair (from R1, R3):
       - Overlapping f(g(y), h(y)) yields the pair (g(y), h(y)).
       - With LPO (h > g), this orients to our first new rule:
         New Rule 1: h(y) -> g(y)

    3. System Simplification and New Pairs:
       - New Rule 1 simplifies R1, creating the pair (f(g(x), g(x)), g(x)).
       - New Rule 1 simplifies R3, creating the pair (f(g(x), g(y)), g(x)).
       - New Rule 1 simplifies the RHS of R2 to g(g(y)).

    4. Processing New Pairs:
       - The pair (f(g(x), g(x)), g(x)) is normalized and oriented.
         - f(g(x), g(x)) -> g(h(g(x)))  (by modified R2)
         - g(h(g(x))) -> g(g(g(x)))  (by New Rule 1)
         - The resulting pair is (g(g(g(x))), g(x)).
         - This orients to our second new rule:
           New Rule 2: g(g(g(x))) -> g(x)
       - The pair (f(g(x), g(y)), g(x)) is oriented. By LPO, f(g(x), g(y)) > g(x) because its
         first argument g(x) is identical to the other term.
         - This gives our third new rule:
           New Rule 3: f(g(x), g(y)) -> g(x)

    5. Termination and Final Ordering:
       - All further critical pairs resolve to trivial ones. The completion is finished.
       - The new rules are:
         - f(g(x), g(y)) -> g(x)
         - g(g(g(x))) -> g(x)
         - h(y) -> g(y)
       - We order them by their left-hand-sides (LHS) using the LPO:
         - f(g(x), g(y)) < g(g(g(x))) < h(y)

    6. Final Result:
       The function prints the final ordered list of rules.
    """
    
    # The final list of rules, ordered by their LHS based on the LPO.
    rule1 = "f(g(x), g(y)) -> g(x)"
    rule2 = "g(g(g(x))) -> g(x)"
    rule3 = "h(y) -> g(y)"
    
    # Print the rules in the correct order, separated by commas.
    print(f"{rule1}, {rule2}, {rule3}")

solve_knuth_bendix()
<<<f(g(x), g(y)) -> g(x), g(g(g(x))) -> g(x), h(y) -> g(y)>>>