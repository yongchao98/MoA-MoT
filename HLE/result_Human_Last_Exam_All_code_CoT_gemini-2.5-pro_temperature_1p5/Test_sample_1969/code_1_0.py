def print_formula_and_numbers():
    """
    This function prints the final derived formula for c_k and explicitly lists the
    numerical constants appearing in it.
    """
    
    # The final derived formula is:
    # c_k = (k * n^2 * (n+1)^2 / 24) * ((3n+1)(n+2) - 6k)
    
    print("The final equation for c_k, the amount for person p_k, is:")
    print("c_k = (k * n**2 * (n + 1)**2 / 24) * ((3*n + 1)*(n + 2) - 6*k)")
    print("\nPrinting each number in the final equation as requested:")
    
    # The numbers are extracted from the formula left to right.
    # Exponent for n
    print(2)
    # Constant in (n+1)
    print(1)
    # Exponent for (n+1)
    print(2)
    # Denominator of the fraction
    print(24)
    # Coefficient of n in (3n+1)
    print(3)
    # Constant in (3n+1)
    print(1)
    # Constant in (n+2)
    print(2)
    # Coefficient of k in -6k (absolute value)
    print(6)

print_formula_and_numbers()