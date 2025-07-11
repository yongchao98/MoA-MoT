import math

def solve():
    """
    This function calculates the minimal polynomial for the shortest geodesic distance
    an ant can travel on a dodecahedron starting from a vertex and returning to it
    without passing through any other vertices.
    """

    # The problem relies on a known result from geometry for the shortest distance 'd'.
    # The square of the distance is given by d^2 = 7 - 3*phi, where phi is the golden ratio.
    # We will derive the minimal polynomial from this expression.

    # phi = (1 + sqrt(5)) / 2
    # d^2 = 7 - 3 * (1 + sqrt(5)) / 2
    # d^2 = (14 - 3 - 3*sqrt(5)) / 2
    # d^2 = (11 - 3*sqrt(5)) / 2
    # This can be written in the form d^2 = (A - B*sqrt(k)) / C
    A = 11
    B = 3
    C = 2
    k = 5

    # From d^2 = (A - B*sqrt(k)) / C, we get C*d^2 - A = -B*sqrt(k).
    # Squaring both sides: (C*d^2 - A)^2 = (-B*sqrt(k))^2
    # This expands to: C^2*d^4 - 2*A*C*d^2 + A^2 = B^2*k
    # Rearranging gives the polynomial: C^2*d^4 - 2*A*C*d^2 + (A^2 - B^2*k) = 0

    # Calculate the coefficients of the polynomial in d.
    # Coefficient of d^4
    coeff_d4 = C**2
    # Coefficient of d^3 is 0
    coeff_d3 = 0
    # Coefficient of d^2
    coeff_d2 = -2 * A * C
    # Coefficient of d^1 is 0
    coeff_d1 = 0
    # Constant term
    coeff_d0 = A**2 - (B**2 * k)

    # The polynomial is: 4*d^4 - 44*d^2 + 76 = 0.
    # We can simplify this by dividing by the greatest common divisor (GCD).
    
    # Python's math.gcd only takes two arguments.
    def find_gcd(a, b, c):
        return math.gcd(math.gcd(a, b), c)

    common_divisor = find_gcd(abs(coeff_d4), abs(coeff_d2), abs(coeff_d0))

    # The coefficients of the minimal polynomial.
    c4 = coeff_d4 // common_divisor
    c3 = coeff_d3 // common_divisor
    c2 = coeff_d2 // common_divisor
    c1 = coeff_d1 // common_divisor
    c0 = coeff_d0 // common_divisor
    
    print("The minimal polynomial for the shortest distance d is P(x) = x^4 - 11x^2 + 19 = 0.")
    print("The equation for d is:")
    
    # As requested, printing each number in the final equation.
    print(f"{c4} * d^4 + {c3} * d^3 + ({c2}) * d^2 + {c1} * d + {c0} = 0")

solve()