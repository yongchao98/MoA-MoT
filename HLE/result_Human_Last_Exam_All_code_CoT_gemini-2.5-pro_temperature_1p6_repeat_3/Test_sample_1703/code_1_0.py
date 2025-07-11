import sympy

def solve_triangle_area():
    """
    Calculates the area of the triangle T(t) as a function of time t.
    """
    # Define t as a symbolic variable
    t = sympy.Symbol('t')
    
    # Radius of the circle and side length of the hexagon
    R = 10
    
    # The area is independent of the hexagon's rotation. We can work in a fixed frame.
    # Vertices of the hexagon. We only need the ones for the sides our triangle is on.
    # The alternating sides are BC, DE, FA.
    B = sympy.Point(R * sympy.cos(sympy.pi/3), R * sympy.sin(sympy.pi/3))
    C = sympy.Point(R * sympy.cos(2*sympy.pi/3), R * sympy.sin(2*sympy.pi/3))
    D = sympy.Point(R * sympy.cos(4*sympy.pi/3), R * sympy.sin(4*sympy.pi/3)) # Note: This is vertex E
    # Wait, the prompt says AB, CD, EF. Let's use those instead for clarity.
    # A = (10, 0), B=(5, 5*sqrt(3))
    # C = (-5, 5*sqrt(3)), D=(-10, 0)
    # E = (-5, -5*sqrt(3)), F=(5, -5*sqrt(3))
    
    # The vertices of T(t) are on sides AB, CD, and EF.
    # Speed v=1. Hexagon side length is R=10.
    # At t=0, the vertices are midpoints (distance 5 from either end of the side).
    # At time t, they have moved a distance of v*t = t from the midpoint.
    # Let's assume they move towards the second vertex of the side (e.g., A->B, C->D, E->F).
    # The distance from the first vertex (A, C, E) is 5 + t.
    
    # Vertex P1(t) on side AB: A=(10,0), B=(5, 5*sqrt(3))
    # Vector B-A = (-5, 5*sqrt(3))
    P1 = sympy.Point(10, 0) + (5 + t)/10 * sympy.Point(-5, 5 * sympy.sqrt(3))
    
    # Vertex P2(t) on side CD: C=(-5, 5*sqrt(3)), D=(-10, 0)
    # Vector D-C = (-5, -5*sqrt(3))
    P2 = sympy.Point(-5, 5 * sympy.sqrt(3)) + (5 + t)/10 * sympy.Point(-5, -5 * sympy.sqrt(3))

    # Calculate the squared distance between P1 and P2 to find the triangle's side length squared.
    # Due to symmetry, the triangle remains equilateral.
    L_squared = P1.distance(P2)**2
    
    # Simplify the expression for L_squared
    L_squared_simple = sympy.simplify(L_squared)
    
    # The area of an equilateral triangle is (sqrt(3)/4) * L^2
    area = (sympy.sqrt(3) / 4) * L_squared_simple
    
    # Extract the coefficients to display the equation clearly.
    # The expression for the squared side length is a quadratic polynomial in t.
    poly_L_sq = sympy.Poly(L_squared_simple, t)
    coeffs = poly_L_sq.all_coeffs() # Should be [a, b, c] for a*t^2 + b*t + c
    
    c_term = coeffs[2] # a*t^2 + b*t + c, so we take the last coeff for constant term
    t_sq_term_coeff = coeffs[0]

    # Print the final result
    print("The area of the triangle T(t) as a function of time t is given by the formula:")
    print(f"Area(t) = (sqrt(3)/4) * ({int(t_sq_term_coeff)}*t**2 + {int(c_term)})")
    
solve_triangle_area()