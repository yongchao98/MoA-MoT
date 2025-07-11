def solve_propositional_logic_problem():
    """
    This script explains the reasoning to find the minimum number of variables
    in a formula psi, which is logically equivalent to a given formula phi.
    The answer is derived symbolically based on the properties of phi.
    """

    # The problem provides a formula φ with n atomic variables (n >= 2).
    # We need to find the minimum number of variables for an equivalent formula ψ.
    # This minimum number is the number of "essential" variables in φ, let's call it m.

    # The key conditions are:
    # 1. φ has 2^(n-1) "discernible" truth-value assignments.
    #    This implies that the number of satisfying assignments for the m essential
    #    variables is k = 2^(n-1).
    # 2. φ is not a tautology.

    # From k = 2^(n-1) and the fact that k must be <= 2^m, we deduce:
    # n - 1 <= m
    # Since m must also be <= n, we have: n - 1 <= m <= n.

    # The case m = n - 1 is ruled out:
    # If m = n-1, the formula must be a tautology on its essential variables to satisfy
    # condition 1. This would make the full formula φ a tautology, violating condition 2.

    # The only possibility left is m = n.
    # This means φ must depend on all n of its variables.
    # Any equivalent formula ψ must therefore also depend on all n variables.
    
    # The final answer is the number of essential variables, which is n.
    final_answer = "n"

    print("The minimum number of distinct atomic variables required in the logically equivalent formula ψ is:")
    print(final_answer)

solve_propositional_logic_problem()