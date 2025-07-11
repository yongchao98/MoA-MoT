def calculate_Dg(g):
    """
    Calculates the value of D_g based on the formula derived by R. Hain.
    D_g is the smallest degree of a finite etale cover of the moduli stack A_g
    of principally polarised abelian varieties of dimension g, such that the
    universal polarisation is represented by a symmetric line bundle on the
    pullback of the universal family.
    """
    if g <= 2:
        # For g=1 and g=2, the degree is 2.
        return 2
    elif g % 2 == 1:  # g is odd and g >= 3
        # For odd g >= 3, the degree is 2^(g-1).
        return 2**(g - 1)
    else:  # g is even and g >= 4
        # For even g >= 4, the degree is 2^g.
        return 2**g

def main():
    """
    Main function to calculate and print the first 4 terms of the sequence D_g.
    """
    # We need the first 4 terms of the sequence, starting from g=1.
    g_values = range(1, 5)
    sequence = []

    print("Calculating the first 4 terms of the sequence D_g using Hain's formula:")
    
    for g in g_values:
        dg = calculate_Dg(g)
        sequence.append(dg)
        if g == 1:
            print(f"For g = {g}, the case is g <= 2. So, D_{g} = 2.")
        elif g == 2:
            print(f"For g = {g}, the case is g <= 2. So, D_{g} = 2.")
        elif g == 3:
            print(f"For g = {g}, g is odd and >= 3. So, D_{g} = 2^({g}-1) = {dg}.")
        elif g == 4:
            print(f"For g = {g}, g is even and >= 4. So, D_{g} = 2^{g} = {dg}.")

    # The final sequence is printed as a string.
    sequence_str = ", ".join(map(str, sequence))
    print(f"\nThe sequence of the first 4 terms of D_g is: {sequence_str}")

if __name__ == "__main__":
    main()