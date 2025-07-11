def solve():
    """
    This function solves the problem by providing the reasoning for each statement
    and then printing the final answer string.
    """
    
    # Reasoning for each statement, leading to the FFFFTT conclusion.
    # We analyze the statements based on two primary counter-example families:
    # 1. Deeply nested strings, like (((...))), where L(x) is proportional to D(x).
    # 2. Widely nested strings, like (()()()...()), where for the outer pair, L(x) is large
    #    while D(x) is small, and for the many inner pairs, L(x) and D(x) are small constants.

    # Statement 1: sum(log L) = O(sum(log D))
    # Test with wide string (()()...()): sum(log L) is Theta(n), but sum(log D) is O(1).
    # False.

    # Statement 2: sum(loglog L) = O(sum(loglog D))
    # Test with wide string: sum(loglog L) is Theta(n), sum(loglog D) is O(1).
    # False.

    # Statement 3: sum(log^5 L) = O(sum(log^5 D))
    # Test with wide string: sum(log^5 L) is Theta(n), sum(log^5 D) is O(1).
    # False.

    # Statement 4: sum(2^sqrt(log L)) = sum(2^O(sqrt(log D)))
    # The RHS notation implies a very strong condition. For the wide string,
    # the LHS grows with n, while the RHS must be a constant, which is a contradiction.
    # False.

    # Statement 5: sum(L^0.1) = O(sum(D^0.11))
    # The exponent on D is larger than on L.
    # For deep strings, L ~ D, so L^0.1 = O(D^0.11).
    # For wide strings, the sums are dominated by the numerous simple () pairs,
    # for which L=2, D=1. The statement holds as n*2^0.1 = O(n*1^0.11).
    # True.

    # Statement 6: sum(L^0.25) = O(sum(D^0.5))
    # Similar to statement 5, the exponent on D is larger than on L (0.5 > 0.25).
    # The same logic applies.
    # True.

    result = "FFFFTT"
    print(result)

solve()
# The final string is a concatenation of the boolean answers (T/F) for the six statements.
# <<<FFFFTT>>>