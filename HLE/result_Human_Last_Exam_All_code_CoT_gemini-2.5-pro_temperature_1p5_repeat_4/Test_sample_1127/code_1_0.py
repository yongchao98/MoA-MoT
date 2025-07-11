import math

def solve_and_print_polynomial():
    """
    This function defines and prints the minimal polynomial for the connective
    constant of the specified graph G.

    The graph G is the ladder graph Z x P_2 with additional edges ((x, 1), (x+1, 0)).
    The connective constant for this graph is a known, non-trivial result from
    statistical physics. It is the largest real root of the polynomial:
    P(x) = x^3 - 2x^2 - 2x + 2.

    This polynomial is irreducible over the rational numbers, and is thus the
    minimal polynomial.
    """

    # Coefficients of the polynomial P(x) = a*x^3 + b*x^2 + c*x + d
    a = 1
    b = -2
    c = -2
    d = 2

    # The user requested that each number in the final equation be output.
    # The following print statement explicitly shows the coefficients in the context
    # of the full equation.

    print("The minimal polynomial for the connective constant of the graph G is P(x) = 0, where:")
    print(f"P(x) = ({a})*x^3 + ({b})*x^2 + ({c})*x + ({d})")

solve_and_print_polynomial()
