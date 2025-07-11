def calculate_minimal_discriminant():
    """
    Calculates the discriminant of the elliptic curve y^2 + xy = x^3 - x - 1.

    This curve is the minimal model for the curve defined by
    y^2 = x^6 + 2*x^3 + 4*x^2 + 4*x + 1.
    The discriminant is calculated using the standard formulas for a general
    Weierstrass equation: y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
    """
    # Coefficients of the minimal model y^2 + xy = x^3 - x - 1
    a1 = 1
    a2 = 0
    a3 = 0
    a4 = -1
    a6 = -1

    # Calculate intermediate quantities b2, b4, b6, b8
    b2 = a1**2 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3**2 + 4 * a6
    b8 = a1**2 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3**2 - a4**2

    # Calculate the terms of the discriminant formula
    # Delta = -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6
    term1 = -b2**2 * b8
    term2 = -8 * b4**3
    term3 = -27 * b6**2
    term4 = 9 * b2 * b4 * b6
    
    # Calculate the final discriminant
    delta = term1 + term2 + term3 + term4

    # Print the explanation and the calculation steps
    print("The minimal Weierstrass model of the curve is y^2 + xy = x^3 - x - 1.")
    print("The coefficients are: a1=1, a2=0, a3=0, a4=-1, a6=-1.")
    print("\nIntermediate values:")
    print(f"b2 = {a1}^2 + 4*({a2}) = {b2}")
    print(f"b4 = 2*({a4}) + {a1}*{a3} = {b4}")
    print(f"b6 = {a3}^2 + 4*({a6}) = {b6}")
    print(f"b8 = {a1}^2*({a6}) + 4*({a2})*({a6}) - {a1}*{a3}*({a4}) + {a2}*{a3}^2 - ({a4})^2 = {b8}")
    
    print("\nThe discriminant Delta is calculated as: Delta = -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6")
    print(f"The equation with the calculated values is:")
    print(f"Delta = -({b2})^2*({b8}) - 8*({b4})^3 - 27*({b6})^2 + 9*({b2})*({b4})*({b6})")
    print(f"Delta = {term1} + {term2} + ({term3}) + {term4}")
    
    print(f"\nThe minimal discriminant is: {delta}")

# Execute the function
calculate_minimal_discriminant()