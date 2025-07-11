def solve_annihilation_time():
    """
    This function prints the formula for the expected time tau.
    The formula is derived by decomposing the problem into contributions
    from three-particle subsystems.
    """
    
    # The formula for the expected time E[tau] is:
    # E[tau] = N1*M1 + M1*N2 + N2*M2 + 0.5 * (N1+M1)*(N2+M2)
    # We will print this formula with explicit coefficients as requested.
    
    # Coefficients for each term
    c1 = 1.0  # for N1*M1
    c2 = 1.0  # for M1*N2
    c3 = 1.0  # for N2*M2
    c4 = 0.5  # for (N1+M1)*(N2+M2)

    formula = (
        f"E[tau] = {c1}*N1*M1 + {c2}*M1*N2 + {c3}*N2*M2 "
        f"+ {c4}*(N1 + M1)*(N2 + M2)"
    )
    
    print(formula)

solve_annihilation_time()