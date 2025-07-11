import math

def solve():
    """
    Computes the simplicial volume of Sigma_31 x Sigma_17.
    """
    g1 = 31
    g2 = 17
    
    # Dimensions of the surfaces
    n1 = 2
    n2 = 2
    
    # Step 1: Compute the simplicial volume of Sigma_g1
    # For a surface of genus g > 1, ||Sigma_g|| = 4g - 4.
    if g1 > 1:
        sv1 = 4 * g1 - 4
        print(f"The simplicial volume of the first surface, Σ_{g1}, is calculated as 4 * {g1} - 4 = {sv1}")
    else:
        sv1 = 0
        print(f"The simplicial volume of the first surface, Σ_{g1}, is {sv1}")

    # Step 2: Compute the simplicial volume of Sigma_g2
    if g2 > 1:
        sv2 = 4 * g2 - 4
        print(f"The simplicial volume of the second surface, Σ_{g2}, is calculated as 4 * {g2} - 4 = {sv2}")
    else:
        sv2 = 0
        print(f"The simplicial volume of the second surface, Σ_{g2}, is {sv2}")
        
    # Step 3: Compute the binomial coefficient for the product formula
    # For a product M_n x M_m, the coefficient is C(n+m, n).
    try:
        binom_coeff = math.comb(n1 + n2, n1)
        print(f"The binomial coefficient for the product, C({n1}+{n2}, {n1}), is {binom_coeff}")
    except ValueError:
        print("Invalid input for binomial coefficient calculation.")
        return

    # Step 4: Compute the final simplicial volume of the product
    # ||Sigma_g1 x Sigma_g2|| = C(n1+n2, n1) * ||Sigma_g1|| * ||Sigma_g2||
    final_sv = binom_coeff * sv1 * sv2
    
    print("\nThe simplicial volume of Σ_{g1} × Σ_{g2} is given by the product of these values.")
    print(f"Final calculation: {binom_coeff} * {sv1} * {sv2} = {final_sv}")
    
    print(f"\nThe simplicial volume of Σ_{g1} × Σ_{g2} is {final_sv}.")

solve()