def calculate_discriminant():
    """
    Calculates the discriminant of the elliptic curve y^2 + y = x^3 - x^2 - 6x - 4.

    This curve is the minimal model of the hyperelliptic curve
    y^2 = x^6 + 2x^3 + 4x^2 + 4x + 1.
    """
    # Coefficients of the general Weierstrass equation
    # y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
    # For y^2 + y = x^3 - x^2 - 6x - 4:
    a1 = 0
    a2 = -1
    a3 = 1
    a4 = -6
    a6 = -4

    # Print the coefficients for clarity
    print("The elliptic curve is y^2 + y = x^3 - x^2 - 6x - 4")
    print("Coefficients:")
    print(f"a1 = {a1}, a2 = {a2}, a3 = {a3}, a4 = {a4}, a6 = {a6}\n")

    # Calculate b-invariants
    b2 = a1**2 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3**2 + 4 * a6
    b8 = a1**2 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3**2 - a4**2

    print("Intermediate b-invariants:")
    print(f"b2 = {b2}")
    print(f"b4 = {b4}")
    print(f"b6 = {b6}")
    print(f"b8 = {b8}\n")

    # Discriminant formula: Delta = -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6
    term1 = -b2**2 * b8
    term2 = -8 * b4**3
    term3 = -27 * b6**2
    term4 = 9 * b2 * b4 * b6
    
    discriminant = term1 + term2 + term3 + term4

    print("Calculating the discriminant Delta = -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6:")
    print(f"Term 1: -({b2})^2 * ({b8}) = {term1}")
    print(f"Term 2: -8 * ({b4})^3 = {term2}")
    print(f"Term 3: -27 * ({b6})^2 = {term3}")
    print(f"Term 4: 9 * ({b2}) * ({b4}) * ({b6}) = {term4}\n")

    print(f"The minimal discriminant is {term1} + {term2} + {term3} + ({term4}) = {discriminant}")

if __name__ == '__main__':
    calculate_discriminant()