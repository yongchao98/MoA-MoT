import sys

def calculate_discriminant():
    """
    Calculates the discriminant of the elliptic curve y^2 + y = x^3 - x^2 - x.
    The general form is y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6.
    """
    
    # Coefficients for the curve y^2 + y = x^3 - x^2 - x
    a1 = 0
    a2 = -1
    a3 = 1
    a4 = -1
    a6 = 0
    
    print("The given curve is birationally equivalent to the elliptic curve y^2 + y = x^3 - x^2 - x.")
    print("This equation is in the general Weierstrass form: y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6")
    print(f"The coefficients are: a1 = {a1}, a2 = {a2}, a3 = {a3}, a4 = {a4}, a6 = {a6}\n")
    
    # Calculate auxiliary quantities b2, b4, b6, b8
    b2 = a1**2 + 4*a2
    b4 = 2*a4 + a1*a3
    b6 = a3**2 + 4*a6
    b8 = a1**2*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3**2 - a4**2

    print("We first compute the auxiliary quantities b2, b4, b6, b8:")
    print(f"b2 = a1^2 + 4*a2 = {a1}^2 + 4*({a2}) = {b2}")
    print(f"b4 = 2*a4 + a1*a3 = 2*({a4}) + {a1}*{a3} = {b4}")
    print(f"b6 = a3^2 + 4*a6 = {a3}^2 + 4*{a6} = {b6}")
    print(f"b8 = a1^2*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3^2 - a4^2 = {a1}^2*{a6} + 4*({a2})*{a6} - {a1}*{a3}*({a4}) + ({a2})*{a3}^2 - ({a4})^2 = {b8}\n")

    # Calculate the discriminant
    delta = -b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
    
    print("The discriminant Delta is calculated using the formula:")
    print("Delta = -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6")
    print(f"Delta = -({b2})^2*({b8}) - 8*({b4})^3 - 27*({b6})^2 + 9*({b2})*({b4})*({b6})")
    
    term1 = -b2**2 * b8
    term2 = -8 * b4**3
    term3 = -27 * b6**2
    term4 = 9 * b2 * b4 * b6
    
    print(f"Delta = -{b2**2}*({b8}) - 8*({b4**3}) - 27*({b6**2}) + {term4}")
    print(f"Delta = {term1} + ({term2}) + ({term3}) + {term4}")
    print(f"Delta = {term1 + term2 + term3 + term4}\n")

    print(f"The minimal discriminant of the curve is {delta}.")

if __name__ == '__main__':
    calculate_discriminant()
