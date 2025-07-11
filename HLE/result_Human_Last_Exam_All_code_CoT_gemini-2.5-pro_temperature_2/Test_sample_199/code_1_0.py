import math

def find_minimal_polynomial_for_distance():
    """
    Calculates the coefficients of the minimal polynomial for the shortest
    geodesic distance on a dodecahedron starting and ending at a vertex.

    The square of the distance, d^2, is known to be 29 - 9*sqrt(5).
    Let y = d^2. The minimal polynomial for y is of the form:
    y^2 - 2*a*y + (a^2 - b^2*k) = 0
    for a value y = a + b*sqrt(k).
    Here, a = 29, b = -9, k = 5.

    Substituting y = d^2 yields a polynomial in d:
    d^4 - 2*a*d^2 + (a^2 - b^2*k) = 0
    """

    # Parameters from d^2 = a + b*sqrt(k)
    a = 29
    b = -9
    k = 5

    # Coefficients for the polynomial for d:
    # c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0 = 0
    c4 = 1
    c3 = 0
    c2 = -2 * a
    c1 = 0
    c0 = a**2 - (b**2) * k

    print("The minimal polynomial for the distance d is of the form:")
    print("x^4 + c2*x^2 + c0 = 0")
    print("\nCalculated coefficients:")
    print(f"Coefficient of x^4: {c4}")
    print(f"Coefficient of x^2: {c2}")
    print(f"Constant term: {c0}")

    # Output the final equation, showing each term
    print("\nThe minimal polynomial equation is:")
    print(f"{c4}*x^4 + {c3}*x^3 + {c2}*x^2 + {c1}*x + {c0} = 0")
    # Simplify for printing
    print("\nSimplified form:")
    print(f"x^4 - {-c2}x^2 + {c0} = 0")


if __name__ == '__main__':
    find_minimal_polynomial_for_distance()
