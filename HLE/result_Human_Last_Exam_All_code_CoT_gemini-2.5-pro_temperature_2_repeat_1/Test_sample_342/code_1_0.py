def solve_sequence():
    """
    This function calculates the first 4 terms of the sequence D_g.
    
    Let g be an integer, and let A_g denote the fine moduli stack over the
    complex numbers of principally polarised abelian varieties (PPAVs) of
    dimension g. Let U_g -> A_g be the universal PPAV.
    
    D_g is the smallest degree of a finite etale cover of A_g such that
    the polarisation is represented by a symmetric line bundle on the
    base-change of U_g to the cover.
    """
    
    # Plan:
    # For g=1, the moduli stack A_1 contains points with non-trivial stabilizers
    # whose action on the 2-torsion points E[2] forces the universal
    # obstruction (Mumford invariant) to be zero. Thus D_1 = 1.
    # For g >= 2, the locus of PPAVs admitting a symmetric representative
    # for the polarization is a proper subvariety of A_g. This implies that
    # the universal obstruction is not identically zero. A cover of degree 2
    # is required to trivialize this obstruction. Thus D_g = 2 for g >= 2.
    
    D = []
    for g in range(1, 5):
        if g == 1:
            d_g = 1
        else:
            d_g = 2
        D.append(d_g)
        
    print(f"The sequence D_g for g=1, 2, 3, 4 is: {', '.join(map(str, D))}")

solve_sequence()