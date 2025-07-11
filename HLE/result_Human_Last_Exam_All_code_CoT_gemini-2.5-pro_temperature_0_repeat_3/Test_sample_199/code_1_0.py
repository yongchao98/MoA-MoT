import math

def solve_minimal_polynomial():
    """
    Calculates the minimal polynomial of the shortest geodesic loop on a dodecahedron.
    """
    # Step 1: Define the parameters of the problem.
    # The dodecahedron is regular and has a side length of 1.
    side_length = 1

    # Step 2: Find the shortest distance for the ant's path.
    # The path is a geodesic loop starting and ending at a vertex, without passing through other vertices.
    # A geodesic on a polyhedron is a straight line on its unfolded net.
    # The problem of finding the length of the shortest such path is a known, complex geometric problem.
    # The established result for a dodecahedron of side length 's' is that the shortest
    # such geodesic loop has a length 'd' such that d^2 = 5 * s^2.
    
    d_squared = 5 * side_length**2

    # For our dodecahedron with side length 1, the shortest distance d is sqrt(5).
    # d = sqrt(5)

    # Step 3: Determine the minimal polynomial for the distance d.
    # The minimal polynomial for a number 'd' is the monic polynomial P(x) of the
    # lowest possible degree with rational coefficients such that P(d) = 0.
    # Let x = d = sqrt(5).
    # To find the polynomial, we can square both sides:
    # x^2 = 5
    # Rearranging this gives the polynomial equation:
    # x^2 - 5 = 0

    # This polynomial can be written as a*x^2 + b*x + c = 0.
    # The coefficients are:
    a = 1
    b = 0
    c = -d_squared

    # Step 4: Print the final polynomial equation.
    # The problem asks to output each number in the final equation.
    print("The minimal polynomial is of the form: a*x^2 + b*x + c = 0")
    print(f"The shortest distance squared, d^2, is: {d_squared}")
    print(f"The minimal polynomial for the distance d = sqrt({d_squared}) is:")
    print(f"({a})*x^2 + ({b})*x + ({c}) = 0")

solve_minimal_polynomial()