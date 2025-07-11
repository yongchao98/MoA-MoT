import math

def solve_minimal_polynomial():
    """
    This function calculates the minimal polynomial of the shortest geodesic loop on a dodecahedron.
    The shortest distance d is known to be sqrt(3) * phi, where phi is the golden ratio.
    """

    # We need to find the minimal polynomial for d.
    # Let's first find the polynomial for x = d^2.
    # x = d^2 = 3 * phi^2 = 3 * (phi + 1)
    # phi = (1 + sqrt(5))/2
    # x = 3 * ( (1+sqrt(5))/2 + 1 ) = 3 * (3+sqrt(5))/2 = (9 + 3*sqrt(5))/2
    # Its conjugate is x_conj = (9 - 3*sqrt(5))/2
    
    # The minimal polynomial for x is z^2 - (x+x_conj)z + (x*x_conj) = 0
    # Coeff for z: -(x+x_conj) = -9
    # Constant term: x*x_conj = (81 - 9*5)/4 = 36/4 = 9
    
    # So, the polynomial for x=d^2 is z^2 - 9z + 9 = 0.
    
    # Substituting z = y^2 gives the polynomial for d: y^4 - 9y^2 + 9 = 0.
    # The polynomial is of the form: c4*y^4 + c3*y^3 + c2*y^2 + c1*y^1 + c0*y^0 = 0
    c4 = 1
    c3 = 0
    c2 = -9
    c1 = 0
    c0 = 9

    print("The minimal polynomial for the shortest distance d is of the form:")
    print("c4*y^4 + c3*y^3 + c2*y^2 + c1*y^1 + c0*y^0 = 0")
    print("\nCalculated coefficients:")
    print(f"c4 = {c4}")
    print(f"c3 = {c3}")
    print(f"c2 = {c2}")
    print(f"c1 = {c1}")
    print(f"c0 = {c0}")

    print("\nThe final equation is:")
    print(f"{c4} * y^4 + {c3} * y^3 + ({c2}) * y^2 + {c1} * y^1 + {c0} * y^0 = 0")

solve_minimal_polynomial()
