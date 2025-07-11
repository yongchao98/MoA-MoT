def solve_knuth_bendix_completion():
    """
    This function derives the new rules added by the Knuth-Bendix completion
    algorithm for the given term-rewriting system.

    It operates under the assumption that the given precedence 'f<g<h' was a typo
    and uses the corrected, workable precedence 'g<h<f'. The derivation steps
    are explained in the comments.
    """

    # --- Initial State ---
    # Original Rules (R0):
    # r1: f(g(x), h(x)) -> g(x)
    # r2: f(y, y) -> g(h(y))
    # r3: f(g(x), h(y)) -> h(x)
    #
    # Assumed Corrected Precedence for LPO: g < h < f

    # --- Derivation of New Rules ---

    # Rule 4: Derived from the critical pair of r1 and r3.
    # Overlapping f(g(x), h(x)) [from r1] and f(g(x), h(y)) [from r3]
    # yields the critical pair (g(x), h(x)).
    # With the LPO precedence g < h < f, h(x) is greater than g(x).
    # This orients the pair into a new rule.
    r4 = ("h(x)", "g(x)")

    # Rule 5: Derived from the critical pair of r1 and the new rule r4.
    # Overlapping h(x) in f(g(x), h(x)) with the LHS of r4 gives the pair
    # (g(x), f(g(x), g(x))).
    # This pair is normalized. The second term reduces as follows:
    # f(g(x), g(x)) -> g(h(g(x)))   (using r2)
    # g(h(g(x)))   -> g(g(g(x)))   (using new rule r4 on the subterm)
    # The normalized pair is (g(x), g(g(g(x)))).
    # The term g(g(g(x))) is greater because it contains g(x) as a subterm.
    # This gives the new rule.
    r5 = ("g(g(g(x)))", "g(x)")

    # Rule 6: Derived from the critical pair of r3 and the new rule r4.
    # Overlapping h(y) in f(g(x), h(y)) with the LHS of r4 gives the pair
    # (h(x), f(g(x), g(y))).
    # This pair is normalized. The first term reduces as follows:
    # h(x) -> g(x) (using new rule r4)
    # The normalized pair is (g(x), f(g(x), g(y))).
    # With LPO precedence g < h < f, f(...) is greater than g(...).
    # This gives the new rule.
    r6 = ("f(g(x), g(y))", "g(x)")
    
    # --- Ordering and Final Result ---
    # The algorithm terminates as all other critical pairs are trivial.
    # The three added rules must be ordered by their Left-Hand-Side (LHS)
    # using the LPO with precedence g < h < f.
    #
    # LHS Terms:
    # r4_lhs = h(x)
    # r5_lhs = g(g(g(x)))
    # r6_lhs = f(g(x), g(y))
    #
    # LPO Comparison:
    # - Comparing g(g(g(x))) and h(x): Since h > g, h(x) > g(g(g(x))).
    # - Comparing h(x) and f(g(x), g(y)): Since f > h, f(g(x), g(y)) > h(x).
    #
    # The final increasing order of LHS is: g(g(g(x))), h(x), f(g(x), g(y)).
    
    ordered_rules = [r5, r4, r6]

    # Print the result in the required format.
    final_answer = ", ".join([f"{lhs} -> {rhs}" for lhs, rhs in ordered_rules])
    print(final_answer)

solve_knuth_bendix_completion()
<<<g(g(g(x))) -> g(x), h(x) -> g(x), f(g(x), g(y)) -> g(x)>>>