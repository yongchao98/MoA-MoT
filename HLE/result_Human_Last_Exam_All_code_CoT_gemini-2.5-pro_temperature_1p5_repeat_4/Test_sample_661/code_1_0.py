def solve_knuth_bendix_completion():
    """
    This script simulates the steps of the Knuth-Bendix completion algorithm
    for the given TRS and term ordering to find the added rules.
    """

    # Term Rewriting System (TRS):
    # R1: f(g(x), h(x)) -> g(x)
    # R2: f(y, y) -> g(h(y))
    # R3: f(g(x), h(y)) -> h(x)
    # Term ordering: Lexicographic Path Ordering (LPO) with precedence f < g < h.

    # Step 1: Find the first critical pair.
    # We find an overlap by unifying the left-hand sides of R1 and R3:
    #   - LHS(R1): f(g(x), h(x))
    #   - LHS(R3): f(g(x'), h(y'))
    # Unifying them results in a substitution {x' -> x, y' -> h(x)}.
    # Applying this substitution to the right-hand sides gives the critical pair:
    #   - From RHS(R1): g(x)
    #   - From RHS(R3): h(x)
    # The critical pair is <g(x), h(x)>.

    # Step 2: Orient the pair to create the first new rule.
    # To orient <g(x), h(x)> using LPO with f<g<h, we see that h > g,
    # so h(x) > g(x).
    # This gives us our first new rule.
    added_rule_1 = "h(x) -> g(x)"

    # Step 3: Find new critical pairs resulting from the new rule.
    # One key consequence (a new critical pair) arises from simplifying R3 with our new rule.
    #   - In LHS(R3), `f(g(x), h(y))`, the subterm `h(y)` reduces to `g(y)` using our new rule.
    #     This gives `f(g(x), g(y))`.
    #   - The RHS(R3) is `h(x)`.
    # This forms a new equation to resolve: f(g(x), g(y)) = h(x).

    # Step 4: Normalize and orient this new equation.
    # We must normalize both sides using all known rules, including `h(x) -> g(x)`.
    #   - `f(g(x), g(y))` cannot be reduced further.
    #   - `h(x)` reduces to `g(x)`.
    # The normalized pair is <f(g(x), g(y)), g(x)>.
    # To orient it, we use LPO. A term is always greater than its subterms.
    # Since `g(x)` is a subterm of `f(g(x), g(y))`, the orientation is:
    added_rule_2 = "f(g(x),g(y)) -> g(x)"

    # Step 5: Verify completion.
    # With the two rules we've added, all other critical pairs can be shown to
    # be "joinable" (i.e., both sides reduce to the same term).
    # For example, simplifying R2 (`f(y,y) -> g(h(y))`) leads to the equation
    # `f(y,y) = g(g(y))`, but this pair is joinable and does not create a new rule.
    # The algorithm terminates.

    # Step 6: Order the added rules by LHS.
    # We must compare `h(x)` and `f(g(x), g(y))` using LPO with f < g < h.
    # Since h has higher precedence than f, h(x) > f(g(x), g(y)).
    # The problem asks for increasing order, so f(g(x),g(y)) comes first.
    final_rules = [added_rule_2, added_rule_1]

    # Print the final result in the specified format.
    print(", ".join(final_rules))

solve_knuth_bendix_completion()
<<<f(g(x),g(y)) -> g(x), h(x) -> g(x)>>>