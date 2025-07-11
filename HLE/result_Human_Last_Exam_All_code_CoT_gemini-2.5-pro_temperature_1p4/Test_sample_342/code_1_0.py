def calculate_Dg(g):
    """
    Calculates the smallest degree D_g of a finite etale cover of A_g
    such that the universal polarisation is represented by a symmetric line bundle.

    Args:
        g: The dimension of the principally polarised abelian varieties (PPAVs).

    Returns:
        The integer value of D_g.
    """
    if g == 1:
        # For g=1, there is an exceptional value for the index [SL(2,Z) : Gamma_1(1,2)], which is 3.
        return 3
    elif g >= 2:
        # For g >= 2, the index [Sp(2g,Z) : Gamma_g(1,2)] is 2^(2g).
        return 2**(2 * g)
    else:
        # g must be a positive integer.
        raise ValueError("The dimension g must be a positive integer.")

def main():
    """
    Calculates and prints the first 4 terms of the sequence D_g.
    """
    print("The sequence of the first 4 terms of D_g is:")
    sequence = []
    for g in range(1, 5):
        dg = calculate_Dg(g)
        sequence.append(dg)
        print(f"D_{g} = {dg}")
    print(f"The sequence is: {', '.join(map(str, sequence))}")

if __name__ == "__main__":
    main()
