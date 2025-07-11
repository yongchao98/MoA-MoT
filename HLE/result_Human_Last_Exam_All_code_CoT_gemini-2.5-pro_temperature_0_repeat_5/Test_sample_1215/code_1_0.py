def solve_logic_puzzle():
    """
    This function solves the logic puzzle by following the reasoning steps outlined.

    1.  The problem asks for the minimum number of variables in a formula `ψ`
        that is logically equivalent to a formula `φ`. This is the same as
        finding the minimum number of *essential variables* that `φ` can have,
        given the conditions.

    2.  Let k be the number of essential variables.

    3.  Lower Bound for k:
        - If k = 0, the formula is a tautology or a contradiction.
        - A tautology is true for 2^n assignments.
        - A contradiction is true for 0 assignments.
        - The problem states φ is true for 2^(n-1) assignments.
        - For n >= 2, 2^(n-1) is not equal to 2^n or 0.
        - Therefore, k must be greater than 0. So, k >= 1.

    4.  Achievable Minimum for k:
        - We need to check if a formula can exist with k = 1.
        - Consider a formula φ that is logically equivalent to a single variable, p_1.
        - This formula is true if and only if p_1 is true. For n variables,
          this means there are 2^(n-1) satisfying assignments (where p_1 is true
          and the other n-1 variables can be anything).
        - This meets the conditions of the problem.
        - This formula has exactly one essential variable (p_1).
        - So, a formula with k = 1 is possible.

    5.  Conclusion:
        - Since k must be at least 1, and we have found a valid case where k = 1,
          the minimum number of essential variables is 1.
    """
    min_variables = 1
    
    # The final equation is simply the value we deduced.
    # We print the number in the final equation as requested.
    print(min_variables)

solve_logic_puzzle()