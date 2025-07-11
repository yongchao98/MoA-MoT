import sympy

def solve_minimal_polynomial():
    """
    Calculates the minimal polynomial for the shortest closed geodesic on a dodecahedron
    starting and ending at a vertex.
    """
    # Step 1: Define the golden ratio, phi, symbolically.
    # The golden ratio phi = (1 + sqrt(5))/2
    phi = (1 + sympy.sqrt(5)) / 2

    # Step 2: Use the known result for the squared shortest distance, d^2.
    # d^2 = 3*phi + 4
    d_squared = 3 * phi + 4

    # Let sympy simplify the expression for d^2
    d_squared_simplified = sympy.simplify(d_squared)
    # This evaluates to (11 + 3*sqrt(5))/2

    # Step 3: Find the minimal polynomial for d^2.
    # We can do this by hand or use sympy's built-in function.
    # Let y = d^2. We have y = (11 + 3*sqrt(5))/2.
    # 2*y - 11 = 3*sqrt(5)
    # Squaring both sides: (2*y - 11)^2 = 45
    # 4*y^2 - 44*y + 121 = 45
    # 4*y^2 - 44*y + 76 = 0
    # y^2 - 11*y + 19 = 0
    
    y = sympy.Symbol('y')
    poly_for_d_squared = sympy.minimal_polynomial(d_squared_simplified, y)
    
    # Step 4: The minimal polynomial for d is found by substituting y = x^2.
    x = sympy.Symbol('x')
    minimal_poly = poly_for_d_squared.subs(y, x**2)

    # Step 5: Extract coefficients and print the final equation.
    # The sympy.Poly object helps to get all coefficients, including zeros.
    p = sympy.Poly(minimal_poly, x)
    coeffs = p.all_coeffs()
    
    # The problem asks to output each number in the final equation.
    print("The minimal polynomial P(x) for the shortest distance is:")
    # The polynomial is x^4 - 11*x^2 + 19 = 0
    print(f"{coeffs[0]}*x^4 + {coeffs[1]}*x^3 + {coeffs[2]}*x^2 + {coeffs[3]}*x + {coeffs[4]} = 0")

solve_minimal_polynomial()