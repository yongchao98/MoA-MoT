def main():
    """
    This script calculates the first four terms of the sequence D_g.
    D_g is the smallest degree of a finite etale cover of the fine moduli stack A_g
    of principally polarised abelian varieties of dimension g, such that the universal
    polarisation can be represented by a symmetric line bundle on the pullback of the
    universal family.
    
    This degree D_g is equal to the number of even theta characteristics for a PPAV
    of dimension g, which is given by the formula: D_g = 2^(g-1) * (2^g + 1).
    """
    
    sequence = []
    
    print("The sequence of the first 4 terms of D_g is calculated as follows:")
    
    for g in range(1, 5):
        # Calculate the components of the formula
        term1_val = 2**(g - 1)
        term2_val = 2**g + 1
        
        # Calculate the final result for D_g
        Dg = term1_val * term2_val
        sequence.append(Dg)
        
        # Print the detailed equation for the current g
        print(f"\nFor g = {g}:")
        print(f"D_{g} = 2^({g}-1) * (2^{g} + 1)")
        print(f"D_{g} = {term1_val} * ({2**g} + 1)")
        print(f"D_{g} = {term1_val} * {term2_val}")
        print(f"D_{g} = {Dg}")
        
    # Print the final sequence
    print("\nThe first 4 terms of the sequence D_g are:")
    print(", ".join(map(str, sequence)))

if __name__ == "__main__":
    main()
