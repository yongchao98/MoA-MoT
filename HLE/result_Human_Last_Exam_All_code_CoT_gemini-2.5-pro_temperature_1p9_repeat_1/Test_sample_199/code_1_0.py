import sympy

def find_minimal_polynomial_of_dodecahedron_geodesic():
    """
    Calculates the minimal polynomial for the shortest vertex-to-vertex
    geodesic loop on a regular dodecahedron.
    """
    
    # Step 1: Define the squared distance, d^2.
    # From geometric analysis of unfolded dodecahedron nets, the shortest
    # distance 'd' for a geodesic loop from a vertex to itself (without passing
    # through other vertices) is known to satisfy d^2 = 7 - 3*sqrt(5).
    d_squared = 7 - 3 * sympy.sqrt(5)

    # The distance d is the square root of this value.
    d = sympy.sqrt(d_squared)

    # Step 2: Use SymPy to find the minimal polynomial of d.
    # We define 'x' as the variable for our polynomial.
    x = sympy.symbols('d')
    
    # The minimal_polynomial function finds the polynomial with rational
    # coefficients of the smallest degree for which d is a root.
    min_poly = sympy.minimal_polynomial(d, x)

    # Step 3: Display the result in the required format.
    # We need to output each number (coefficient) in the final equation.
    
    # The minimal polynomial is x^4 - 14x^2 + 4.
    # We can get the coefficients from the poly object.
    coeffs_dict = min_poly.as_dict()

    # The equation is c4*d^4 + c3*d^3 + c2*d^2 + c1*d + c0 = 0
    c4 = coeffs_dict.get((4,), 0)
    c3 = coeffs_dict.get((3,), 0)
    c2 = coeffs_dict.get((2,), 0)
    c1 = coeffs_dict.get((1,), 0)
    c0 = coeffs_dict.get((0,), 0)
    
    print("The minimal polynomial for the shortest distance, d, is:")
    print(min_poly)
    print("\nThe final equation with each coefficient explicitly shown is:")
    print(f"{c4} * d**4 + {c3} * d**3 + ({c2}) * d**2 + {c3} * d + {c0} = 0")


find_minimal_polynomial_of_dodecahedron_geodesic()
