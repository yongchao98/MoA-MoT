def solve():
    """
    This function prints the complexity classes for languages G and G'.

    G = {M | exists w forall s: M accepts sw}
    G' = {(k,M) | exists w forall s: |s|=k implies M accepts sw}

    Based on the analysis:
    - Both G and G' are NP-complete.
    - NP is equivalent to Sigma_1^P in the polynomial hierarchy.
    - Therefore, the lowest rung for both is NP.
    """
    # The complexity classes for G and G' respectively.
    # The format is "CLASS_G, CLASS_G_PRIME"
    print("NP, NP")

solve()