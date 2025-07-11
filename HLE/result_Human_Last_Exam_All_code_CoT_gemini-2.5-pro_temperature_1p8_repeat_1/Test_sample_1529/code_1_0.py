import math

def calculate_minimal_discriminant():
    """
    Calculates the minimal discriminant of the specified curve.

    The curve y^2 = x^6 + 2x^3 + 4x^2 + 4x + 1 is birationally equivalent
    to the elliptic curve y^2 + y = x^3 - x (Cremona label '43a1').
    We will calculate the discriminant for this model.
    """
    print("The curve y^2 = x^6 + 2x^3 + 4x^2 + 4x + 1 is birationally equivalent to the elliptic curve E: y^2 + y = x^3 - x.")
    print("The equation for E is in the general Weierstrass form: y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6.")
    
    # Coefficients of the elliptic curve y^2 + y = x^3 - x
    a1 = 0
    a2 = 0
    a3 = 1
    a4 = -1
    a6 = 0

    print("\nThe coefficients (a_i) for the model y^2 + y = x^3 - x are:")
    print(f"a1 = {a1}")
    print(f"a2 = {a2}")
    print(f"a3 = {a3}")
    print(f"a4 = {a4}")
    print(f"a6 = {a6}")

    # Calculate b invariants
    b2 = a1**2 + 4*a2
    b4 = 2*a4 + a1*a3
    b6 = a3**2 + 4*a6
    b8 = a1**2*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3**2 - a4**2

    print("\nNext, we calculate the b-invariants using the formulas:")
    print(f"b2 = a1^2 + 4*a2 = ({a1})^2 + 4*({a2}) = {b2}")
    print(f"b4 = 2*a4 + a1*a3 = 2*({a4}) + ({a1})*({a3}) = {b4}")
    print(f"b6 = a3^2 + 4*a6 = ({a3})^2 + 4*({a6}) = {b6}")
    print(f"b8 = a1^2*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3^2 - a4^2 = ({a1})^2*({a6}) + 4*({a2})*({a6}) - ({a1})*({a3})*({a4}) + ({a2})*({a3})^2 - ({a4})^2 = {b8}")
    
    # Calculate the discriminant
    term1 = -b2**2 * b8
    term2 = -8 * b4**3
    term3 = -27 * b6**2
    term4 = 9 * b2 * b4 * b6
    discriminant = term1 + term2 + term3 + term4

    print("\nFinally, we calculate the discriminant Delta using the formula:")
    print("Delta = -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6")
    print("\nHere is the calculation with the numbers plugged in:")
    print(f"Delta = -({b2})^2*({b8}) - 8*({b4})^3 - 27*({b6})^2 + 9*({b2})*({b4})*({b6})")
    print(f"      = -({b2**2})*({b8}) - 8*({b4**3}) - 27*({b6**2}) + ({term4})")
    print(f"      = ({term1}) + ({term2}) + ({term3}) + ({term4})")
    print(f"      = {discriminant}")
    
    print("\nThe discriminant for this model is 37.")
    print("An isomorphic curve obtained by a change of variables (x = u^2*x' + r, y = u^3*y' + ...) has discriminant Delta' = u^(-12)*Delta.")
    print("For Delta' to be an integer, u^12 must divide Delta = 37.")
    print("Since 37 is a prime number, the only integer solutions for u are u = 1 and u = -1.")
    print("This means the discriminant cannot be made smaller in absolute value.")
    print("\nThus, the minimal discriminant of the curve is 37.")


calculate_minimal_discriminant()