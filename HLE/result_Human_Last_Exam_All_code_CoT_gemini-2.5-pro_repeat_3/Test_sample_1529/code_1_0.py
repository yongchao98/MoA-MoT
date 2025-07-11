import math

def solve_discriminant():
    """
    Calculates the minimal discriminant of the elliptic curve y^2 + y = x^3 - x,
    which is associated with the hyperelliptic curve y^2 = x^6 + 2x^3 + 4x^2 + 4x + 1.
    """
    
    # The elliptic curve is y^2 + y = x^3 - x.
    # This is in the generalized Weierstrass form:
    # y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
    
    # Step 1: Identify coefficients
    a1 = 0
    a2 = 0
    a3 = 1
    a4 = -1
    a6 = 0
    
    # Step 2: Calculate b-invariants
    b2 = a1**2 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3**2 + 4 * a6
    
    # Step 3: Calculate c-invariants
    c4 = b2**2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    
    # Step 4: Calculate the discriminant using the formula
    # Delta = (c4^3 - c6^2) / 1728
    # We use integer division // as the result is guaranteed to be an integer.
    c4_cubed = c4**3
    c6_squared = c6**2
    numerator = c4_cubed - c6_squared
    denominator = 1728
    discriminant = numerator // denominator

    # Output the numbers in the final equation and the result
    print("The minimal discriminant is calculated using the formula: Delta = (c4^3 - c6^2) / 1728")
    print(f"The value for c4 is: {c4}")
    print(f"The value for c6 is: {c6}")
    print(f"The equation is: Delta = ({c4}^3 - ({c6})^2) / {denominator}")
    print(f"Which evaluates to: Delta = ({c4_cubed} - {c6_squared}) / {denominator}")
    print(f"The final result is: {discriminant}")

solve_discriminant()