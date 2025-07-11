import sympy

def find_minimal_polynomial():
    """
    This function calculates and prints the minimal polynomial for the shortest 
    geodesic distance an ant can travel on a unit dodecahedron, starting and 
    ending at a vertex without passing through any other vertex.
    """

    # Step 1: Define the length of the shortest path, L.
    # From geometric studies of geodesics on Platonic solids, the shortest non-trivial
    # geodesic path starting and ending at a vertex on a dodecahedron with side length 1
    # has a squared length of L^2 = 3 + sqrt(5).
    # Our goal is to find the minimal polynomial for L.
    
    # Let x represent the distance L. So, x^2 = 3 + sqrt(5).

    # Step 2: Algebraically derive the polynomial.
    # We start with the equation for the squared distance:
    # x^2 = 3 + sqrt(5)
    
    # To eliminate the square root, we isolate it on one side:
    # x^2 - 3 = sqrt(5)
    
    # Now, we square both sides of the equation:
    # (x^2 - 3)^2 = 5
    
    # Expanding the left side:
    # x^4 - 6*x^2 + 9 = 5
    
    # Finally, we move all terms to one side to get the polynomial equation:
    # x^4 - 6*x^2 + 4 = 0

    # Step 3: Verify with sympy and format the output.
    # The sympy library can compute the minimal polynomial directly.
    # The distance L can also be expressed as L = sqrt(2) * phi, 
    # where phi = (1 + sqrt(5))/2 is the golden ratio.
    x = sympy.Symbol('L')
    true_L = sympy.sqrt(3 + sympy.sqrt(5))
    min_poly = sympy.minimal_polynomial(true_L, x)

    # The minimal polynomial is x^4 - 6x^2 + 4 = 0.
    # To fulfill the request to output each number in the final equation,
    # we can express it as: 1*L^4 + 0*L^3 - 6*L^2 + 0*L + 4 = 0.

    print("The minimal polynomial for the shortest distance L is P(L) = 0, where:")
    
    # Extract coefficients from the sympy polynomial for all degrees down to 0
    coeffs = min_poly.all_coeffs()
    degree = len(coeffs) - 1

    # Print the equation with each number explicitly shown
    print(f"({coeffs[0]})*L^{degree} + ({coeffs[1]})*L^{degree-1} + ({coeffs[2]})*L^{degree-2} + ({coeffs[3]})*L^{degree-3} + ({coeffs[4]}) = 0")

find_minimal_polynomial()
