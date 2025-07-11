def solve_sequence():
    """
    Calculates the first 4 terms of the sequence D_g.
    
    The value D_g is the smallest degree of a finite etale cover of the fine moduli stack A_g 
    of principally polarised abelian varieties (PPAVs) of dimension g, such that the 
    polarisation is represented by a symmetric line bundle on the base-change of the universal 
    PPAV to the cover.

    The value of D_g depends on the parity of g:
    - If g is odd, D_g = 1.
    - If g is even, D_g = 2.
    """
    
    sequence = []
    # We calculate the terms for g = 1, 2, 3, 4.
    for g in range(1, 5):
        if g % 2 == 1:
            # g is odd
            d_g = 1
        else:
            # g is even
            d_g = 2
        sequence.append(d_g)
    
    # The final instruction asks to output each number, which we interpret
    # as printing the final sequence.
    print(f"The first 4 terms of the sequence D_g are: {', '.join(map(str, sequence))}")

solve_sequence()