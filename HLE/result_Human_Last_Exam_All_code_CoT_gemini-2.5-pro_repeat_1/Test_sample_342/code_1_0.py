def solve_sequence():
    """
    Calculates the first 4 terms of the sequence D_g.

    Let g be an integer. A_g is the fine moduli stack of g-dimensional
    principally polarised abelian varieties (PPAVs). U_g is the universal
    PPAV over A_g. D_g is the smallest degree of a finite etale cover of A_g
    such that the universal polarisation on the pullback of U_g is
    represented by a symmetric line bundle.

    The value of D_g depends on the parity of g. According to known results
    in the theory of moduli of abelian varieties:
    - D_g = 1 if g is odd.
    - D_g = 2 if g is even.
    """

    # We need to compute D_g for g = 1, 2, 3, 4.
    sequence = []
    for g in range(1, 5):
        if g % 2 == 1:
            # g is odd
            Dg = 1
        else:
            # g is even
            Dg = 2
        sequence.append(Dg)

    # Print the result in the format D_g = value
    for g, val in enumerate(sequence, 1):
        print(f"D_{g} = {val}")

solve_sequence()
