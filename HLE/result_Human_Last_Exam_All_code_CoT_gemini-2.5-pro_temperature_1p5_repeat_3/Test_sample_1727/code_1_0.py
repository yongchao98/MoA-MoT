import math

def print_chi_formula():
    """
    This function prints the derived symbolic formula for the amplitude factor chi.
    """
    
    print("The stable amplitude factor chi is a function of the dimensionless parameter K = k*h, where k is the wavenumber and h is the plate thickness.")
    print("Let C = cosh(K) and S = sinh(K).")
    
    # The derived formula is chi = Numerator / Denominator
    
    # Numerator of the expression for chi
    numerator_str = "-2 * K * C"
    
    # Denominator parts of the expression for chi
    denominator_part1_str = "K * C * S"
    denominator_part2_str = "2 * K**2"
    denominator_part3_str = "2 * C**2"
    
    print("\nThe final equation for the amplitude factor chi is:")
    print(f"chi = ( {numerator_str} ) / ( {denominator_part1_str} + {denominator_part2_str} + {denominator_part3_str} )")

    print("\nBreaking down the equation into its components:")
    
    print("\nNumerator:")
    print("N = -2 * K * cosh(K)")
    
    print("\nDenominator:")
    print("D = K * cosh(K) * sinh(K) + 2 * K^2 + 2 * cosh(K)^2")

print_chi_formula()