def solve():
    """
    This function determines which of the 12 given identities necessarily follow from the problem's assumption.

    The logic is as follows:
    1.  The assumption Psi(k;l;m) = 0 is analyzed to derive its consequences.
    2.  General identities that hold true based on the algebraic structure (commutative, idempotent monoid M acting on an abelian group G) are identified.
    3.  Each of the 12 statements is evaluated. A statement "follows necessarily" if it is a consequence of the assumption or if it is a general identity.

    -   Statements 1, 2, 3, 4, 5, 9, 10 are found to be not necessarily true by considering counterexamples or showing that they do not follow from the derived consequences.
    -   Statement 6, k.Phi^2(l;m) = 0, is proven to be a direct consequence of the assumption.
    -   Statements 7, (lm).Phi^2(k;m) = 0, and 8, (klm).Phi^2(k;l) = 0, are proven to be general identities stemming from the idempotent property of M.
    -   Statements 11, (lm).Phi^3(k;l;m) = 0, and 12, (klm).Phi^3(k;l;m) = 0, are also proven to be general identities.

    Therefore, the necessarily true identities are 6, 7, 8, 11, and 12.
    """
    # The numbers of the identities that are necessarily true.
    true_identities = [6, 7, 8, 11, 12]

    # The result is formatted as a comma-separated string, in increasing order.
    result = ",".join(map(str, sorted(true_identities)))
    print(result)

solve()