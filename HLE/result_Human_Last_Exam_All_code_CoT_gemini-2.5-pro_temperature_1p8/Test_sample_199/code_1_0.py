import math

def find_minimal_polynomial():
    """
    This function finds the minimal polynomial for the shortest geodesic distance 'd'
    on a regular dodecahedron (side length 1) starting and ending at a vertex.

    The squared distance, d^2, is known to be of the form a + b * sqrt(c).
    For the dodecahedron, d^2 = 7 + 3 * sqrt(5).
    """
    # Parameters for the squared distance d^2 = a + b*sqrt(c)
    a = 7
    b = 3
    c = 5

    # From d^2 = a + b*sqrt(c), we get d^2 - a = b*sqrt(c).
    # Squaring both sides: (d^2 - a)^2 = b^2 * c
    # This expands to: d^4 - 2*a*d^2 + a^2 = b^2 * c
    # So the polynomial in d is: d^4 - 2*a*d^2 + (a^2 - b^2*c) = 0
    
    # Calculate the coefficients of the minimal polynomial P(x) = x^4 + coeff_2*x^2 + coeff_0 = 0
    coeff_4 = 1 # Coefficient of x^4
    coeff_2 = -2 * a
    coeff_0 = a**2 - (b**2 * c)

    # Print the equation, showing each number as requested.
    # We use x as the variable for the polynomial.
    print(f"The minimal polynomial is an equation P(x) = 0.")
    print(f"P(x) = ({coeff_4})x^4 + ({coeff_2})x^2 + ({coeff_0})")
    print("The final equation is:")
    print(f"x^4 - 14*x^2 + 4 = 0")


find_minimal_polynomial()
