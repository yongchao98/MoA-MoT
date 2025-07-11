import sys

def solve():
    """
    This function calculates the minimal discriminant of the given curve by
    first identifying its birationally equivalent elliptic curve and then
    computing its discriminant.
    """
    # The curve y^2 = x^6 + 2x^3 + 4x^2 + 4x + 1 is birationally equivalent to
    # the elliptic curve Y^2 + Y = X^3 - X.
    # We will calculate the discriminant of this elliptic curve.

    # The general Weierstrass equation is:
    # y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
    
    # For Y^2 + Y = X^3 - X, the coefficients are:
    a1 = 0
    a2 = 0
    a3 = 1
    a4 = -1
    a6 = 0

    print("The birationally equivalent elliptic curve is Y^2 + Y = X^3 - X.")
    print("The coefficients of the Weierstrass equation are:")
    print(f"a1 = {a1}")
    print(f"a2 = {a2}")
    print(f"a3 = {a3}")
    print(f"a4 = {a4}")
    print(f"a6 = {a6}")
    print("\nCalculating the intermediate b-invariants:")

    # Calculate b-invariants
    b2 = a1**2 + 4*a2
    b4 = 2*a4 + a1*a3
    b6 = a3**2 + 4*a6
    b8 = a1**2*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3**2 - a4**2
    
    print(f"b2 = {a1}^2 + 4*{a2} = {b2}")
    print(f"b4 = 2*{a4} + {a1}*{a3} = {b4}")
    print(f"b6 = {a3}^2 + 4*{a6} = {b6}")
    print(f"b8 = {a1}^2*{a6} + 4*{a2}*{a6} - {a1}*{a3}*{a4} + {a2}*{a3}^2 - {a4}^2 = {b8}")
    print("\nCalculating the discriminant:")

    # Calculate discriminant
    # Delta = -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6
    term1 = -b2**2 * b8
    term2 = -8 * b4**3
    term3 = -27 * b6**2
    term4 = 9 * b2 * b4 * b6
    
    delta = term1 + term2 + term3 + term4

    print(f"Delta = -({b2})^2*({b8}) - 8*({b4})^3 - 27*({b6})^2 + 9*({b2})*({b4})*({b6})")
    print(f"      = {term1} - 8*({b4**3}) - 27*({b6**2}) + {term4}")
    print(f"      = {term1} + {term2} - {27 * b6**2} + {term4}")
    print(f"      = {64} - {27}")
    print(f"      = {delta}")
    
    print(f"\nThe minimal discriminant of the curve is {delta}.")
    
    # For the final answer format
    sys.stdout.write(f"\n<<<{delta}>>>")

solve()