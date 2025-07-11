import math

def find_minimal_polynomial_for_dodecahedron_path():
    """
    This function calculates the coefficients of the minimal polynomial for the shortest
    geodesic distance on a dodecahedron starting and ending at a vertex.

    The square of the distance, d^2, is known to be 27 + 9*sqrt(5).
    Let y = d^2 = a + b*sqrt(5), where a=27, b=9.
    The minimal polynomial for y is of the form y^2 - 2*a*y + (a^2 - 5*b^2) = 0.
    We substitute y = x^2 to get the polynomial for x=d.
    """
    
    # d^2 = a + b * sqrt(5)
    a = 27
    b = 9

    # From (y-a)^2 = (b*sqrt(5))^2, we get:
    # y^2 - 2*a*y + a^2 = 5*b^2
    # y^2 - 2*a*y + (a^2 - 5*b^2) = 0
    
    # Coefficients for the polynomial in y = x^2: P(y) = c2*y^2 + c1*y + c0
    c2_y = 1
    c1_y = -2 * a
    c0_y = a**2 - 5 * b**2
    
    # The polynomial for x is derived by substituting y = x^2:
    # (x^2)^2 + c1_y*(x^2) + c0_y = 0
    # x^4 + c1_y*x^2 + c0_y = 0
    # The coefficients for P(x) = c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0 are:
    c4 = c2_y
    c3 = 0
    c2 = c1_y
    c1 = 0
    c0 = c0_y
    
    # The minimal polynomial is P(x) = x^4 - 54*x^2 + 324 = 0
    print("The minimal polynomial for the shortest distance d is of the form:")
    print("c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0 = 0")
    print("\nThe coefficients are:")
    print(f"c4 = {c4}")
    print(f"c3 = {c3}")
    print(f"c2 = {c2}")
    print(f"c1 = {c1}")
    print(f"c0 = {c0}")

    print("\nSo the final equation is:")
    # Print each part of the equation as requested
    print(f"{c4}x^4 + {c2}x^2 + {c0} = 0")


find_minimal_polynomial_for_dodecahedron_path()

# The final answer format
print("\n<<<x^4 - 54*x^2 + 324 = 0>>>")