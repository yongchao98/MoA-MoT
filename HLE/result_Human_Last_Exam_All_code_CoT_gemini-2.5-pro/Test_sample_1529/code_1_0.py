def calculate_minimal_discriminant():
    """
    Calculates the minimal discriminant of the curve y^2 = x^6 + 2x^3 + 4x^2 + 4x + 1.

    The curve is birationally equivalent to the elliptic curve y^2 = x^3 - x^2 - x.
    This function computes the discriminant of this elliptic curve.
    """
    
    # Coefficients of the elliptic curve y^2 = x^3 + a2*x^2 + a4*x + a6
    a2 = -1
    a4 = -1
    a6 = 0
    
    print(f"The associated elliptic curve is y^2 = x^3 + ({a2})x^2 + ({a4})x + ({a6})")
    print("Calculating the discriminant using the standard formula.")

    # Intermediate quantities
    b2 = 4 * a2
    b4 = 2 * a4
    b6 = 4 * a6
    b8 = 4 * a2 * a6 - a4**2
    
    print(f"b2 = 4 * a2 = 4 * ({a2}) = {b2}")
    print(f"b4 = 2 * a4 = 2 * ({a4}) = {b4}")
    print(f"b6 = 4 * a6 = 4 * ({a6}) = {b6}")
    print(f"b8 = 4 * a2 * a6 - a4^2 = 4 * ({a2}) * ({a6}) - ({a4})^2 = {b8}")

    # Discriminant formula: Delta = -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6
    term1 = -b2**2 * b8
    term2 = -8 * b4**3
    term3 = -27 * b6**2
    term4 = 9 * b2 * b4 * b6
    
    discriminant = term1 + term2 + term3 + term4
    
    print("\nCalculating the terms of the discriminant formula: Delta = -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6")
    print(f"-b2^2*b8 = -({b2})^2 * ({b8}) = {term1}")
    print(f"-8*b4^3 = -8 * ({b4})^3 = {term2}")
    print(f"-27*b6^2 = -27 * ({b6})^2 = {term3}")
    print(f"9*b2*b4*b6 = 9 * ({b2}) * ({b4}) * ({b6}) = {term4}")

    print(f"\nFinal Discriminant Delta = {term1} + {term2} + {term3} + {term4} = {discriminant}")

calculate_minimal_discriminant()