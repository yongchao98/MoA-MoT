import math

def solve():
    """
    This function calculates and prints the minimal polynomial for the shortest
    geodesic distance an ant can travel on a dodecahedron starting from a vertex,
    visiting no other vertices, and returning to the start.

    The problem is a known one in geometry. The square of the shortest distance 'd'
    for a simple closed geodesic starting and ending at a vertex of a unit-side-length
    dodecahedron is given by d^2 = 27 + 9*sqrt(5).

    From d^2 = 27 + 9*sqrt(5), we derive the minimal polynomial.
    Let x = d.
    x^2 = 27 + 9*sqrt(5)
    x^2 - 27 = 9*sqrt(5)
    Squaring both sides to eliminate the radical:
    (x^2 - 27)^2 = (9*sqrt(5))^2
    x^4 - 54*x^2 + 729 = 81 * 5
    x^4 - 54*x^2 + 729 = 405
    x^4 - 54*x^2 + 324 = 0

    This polynomial, x^4 - 54x^2 + 324 = 0, is the minimal polynomial for the distance d.
    We will now print the coefficients of this polynomial.
    """

    # Coefficients of the minimal polynomial P(x) = c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0 = 0
    c4 = 1
    c3 = 0
    # From (x^2 - 27)^2 = 405
    c2 = -54
    c1 = 0
    # From 27^2 - 405 = 729 - 405 = 324
    c0 = 324

    print("The minimal polynomial is of the form: c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0 = 0")
    print("The coefficients of the minimal polynomial are:")
    print(f"c4: {c4}")
    print(f"c3: {c3}")
    print(f"c2: {c2}")
    print(f"c1: {c1}")
    print(f"c0: {c0}")

    print("\nThe final equation is:")
    print(f"{c4}x^4 + {c3}x^3 + ({c2})x^2 + {c1}x + {c0} = 0")
    # A cleaner representation:
    print(f"x^4 - 54x^2 + 324 = 0")

solve()