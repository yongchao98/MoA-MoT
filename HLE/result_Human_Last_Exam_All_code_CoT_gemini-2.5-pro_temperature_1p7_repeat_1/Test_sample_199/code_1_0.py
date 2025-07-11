import math

def solve():
    """
    This function calculates and prints the minimal polynomial for the shortest 
    geodesic distance an ant can walk from a vertex back to itself on a regular
    dodecahedron without passing through other vertices.

    The squared length of the shortest path is d^2 = 13 + 6*sqrt(5).
    We find the minimal polynomial for d.
    """
    
    # The problem is to find the minimal polynomial of d. Let x = d.
    # x^2 = 13 + 6*sqrt(5)
    # To eliminate the square root, we rearrange and square the equation.
    # (x^2 - 13)^2 = (6*sqrt(5))^2
    # x^4 - 26*x^2 + 169 = 36 * 5
    # x^4 - 26*x^2 + 169 = 180
    # x^4 - 26*x^2 - 11 = 0
    
    # Coefficients of the polynomial p(x) = c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0 = 0
    c4 = 1
    c3 = 0
    c2 = -26
    c1 = 0
    c0 = -11
    
    print("The shortest non-trivial geodesic distance 'd' for an ant starting at a vertex of a regular dodecahedron of side length 1 and returning to it without passing through any other vertex satisfies the equation:")
    print(f"{c4}*d^4 + ({c2})*d^2 + ({c0}) = 0")
    print("So the minimal polynomial is:")
    print(f"x^4 - 26*x^2 - 11 = 0")

solve()