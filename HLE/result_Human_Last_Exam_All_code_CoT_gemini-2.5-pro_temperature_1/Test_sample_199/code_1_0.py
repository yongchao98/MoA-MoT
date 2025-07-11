import math

def solve_minimal_polynomial():
    """
    This function calculates the minimal polynomial for the shortest geodesic distance
    on a dodecahedron starting and ending at the same vertex.

    The problem is based on a known result from geometry. The square of the shortest
    distance, d^2, is given by a value of the form a + b*sqrt(c).
    From this, we derive a polynomial equation for d with integer coefficients.
    """

    # The square of the shortest distance (d^2) is known to be 13 - 3*sqrt(5).
    # Let d^2 = y = a + b * sqrt(c)
    a = 13
    b = -3
    c = 5

    # We have y = a + b*sqrt(c).
    # To get a polynomial with rational coefficients, we isolate the square root:
    # y - a = b*sqrt(c)
    # Then we square both sides:
    # (y - a)^2 = (b*sqrt(c))^2
    # y^2 - 2*a*y + a^2 = b^2*c
    # y^2 - 2*a*y + (a^2 - b^2*c) = 0

    # Let's calculate the coefficients for this polynomial in y.
    coeff_y = -2 * a
    constant_term_y = a**2 - (b**2 * c)

    # Now, we substitute y = d^2 (or x^2 in the polynomial) back.
    # The polynomial for d (let's use 'x' as the variable) is:
    # x^4 + coeff_y * x^2 + constant_term_y = 0
    
    coeff_x2 = int(coeff_y)
    constant_term_x = int(constant_term_y)

    # Print the explanation and the final equation.
    print("The problem is to find the minimal polynomial of the shortest geodesic distance (d) on a dodecahedron.")
    print("The path starts and ends at the same vertex without passing through any other vertices.")
    print(f"From geometry, the square of this distance is known: d^2 = {a} + {b}*sqrt({c}).")
    print("\nTo find the minimal polynomial for d, we let x = d and start with x^2 = 13 - 3*sqrt(5).")
    print("We rearrange to x^2 - 13 = -3*sqrt(5).")
    print("Squaring both sides gives (x^2 - 13)^2 = (-3*sqrt(5))^2, which is x^4 - 26*x^2 + 169 = 45.")
    print("This simplifies to the final minimal polynomial equation:\n")
    
    # We output each number in the final equation as requested.
    # The format is x^4 + (coeff_x2)x^2 + (constant_term_x) = 0
    sign_coeff = '+' if coeff_x2 > 0 else '-'
    sign_const = '+' if constant_term_x > 0 else '-'
    
    print(f"x^4 {sign_coeff} {abs(coeff_x2)}*x^2 {sign_const} {abs(constant_term_x)} = 0")


solve_minimal_polynomial()
