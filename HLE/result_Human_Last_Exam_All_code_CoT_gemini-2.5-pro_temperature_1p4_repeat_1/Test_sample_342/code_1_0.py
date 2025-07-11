def compute_D_g_sequence():
    """
    Computes and prints the first 4 terms of the sequence D_g.

    D_g is the smallest degree of a finite etale cover of the fine moduli stack A_g
    of principally polarised abelian varieties (PPAVs) of dimension g, such that
    the universal polarisation can be represented by a symmetric line bundle.
    """

    def get_Dg(g):
        """
        Calculates the value of D_g based on mathematical theory.

        The ability to represent the universal polarisation by a symmetric line bundle
        is obstructed by a class in the second cohomology group H^2(Sp(2g, Z), Z/2Z).

        Args:
            g (int): The dimension of the abelian varieties.

        Returns:
            int: The value of D_g.
        """
        if g == 1:
            # For g=1, the group is SL(2,Z), and its second cohomology with Z/2Z
            # coefficients is trivial. H^2(SL(2,Z), Z/2Z) = 0.
            # This means there is no obstruction, so no cover is needed.
            # The degree of the cover is 1.
            return 1
        elif g >= 2:
            # For g >= 2, the group is Sp(2g,Z), and its second cohomology group
            # H^2(Sp(2g,Z), Z/2Z) = Z/2Z.
            # The obstruction is non-trivial, and its order is 2. To kill an
            # obstruction of order 2, a double cover is necessary and sufficient.
            # The minimal degree of the cover is 2.
            return 2
        else:
            raise ValueError("Dimension g must be a positive integer.")

    # Calculate the sequence for g = 1, 2, 3, 4
    sequence = [get_Dg(g) for g in range(1, 5)]

    # Print the numbers of the sequence, as requested.
    print(", ".join(map(str, sequence)))

compute_D_g_sequence()