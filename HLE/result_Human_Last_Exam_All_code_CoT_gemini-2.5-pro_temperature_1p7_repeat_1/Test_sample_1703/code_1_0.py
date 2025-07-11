import sympy

def calculate_triangle_area():
    """
    Calculates the area of the triangle T(t) as a function of time t.
    """
    # Define symbols for time and mathematical constants.
    t = sympy.Symbol('t', real=True, nonnegative=True)
    R = 10
    sqrt3 = sympy.sqrt(3)

    # Step 1: Define the hexagon vertices in a convenient orientation.
    # The center is at (0,0), radius R=10.
    # We choose an orientation with two vertices on the y-axis for simplicity.
    V1 = sympy.Point2D(R * sympy.cos(sympy.pi / 2), R * sympy.sin(sympy.pi / 2))     # (0, 10)
    V2 = sympy.Point2D(R * sympy.cos(sympy.pi / 6), R * sympy.sin(sympy.pi / 6))     # (5*sqrt(3), 5)
    V3 = sympy.Point2D(R * sympy.cos(11 * sympy.pi / 6), R * sympy.sin(11 * sympy.pi / 6)) # (5*sqrt(3), -5)
    V4 = sympy.Point2D(R * sympy.cos(3 * sympy.pi / 2), R * sympy.sin(3 * sympy.pi / 2)) # (0, -10)
    V5 = sympy.Point2D(R * sympy.cos(7 * sympy.pi / 6), R * sympy.sin(7 * sympy.pi / 6))  # (-5*sqrt(3), -5)
    V6 = sympy.Point2D(R * sympy.cos(5 * sympy.pi / 6), R * sympy.sin(5 * sympy.pi / 6))  # (-5*sqrt(3), 5)
    
    # Alternating sides for the triangle vertices are V6V1, V2V3, V4V5.
    
    # Step 2: Parameterize the position of the triangle's vertices P, Q, R.
    # Vertices start at midpoints (t=0) and move with speed v=1.
    # Hexagon side length L=R=10.
    # The distance traveled from the midpoint along a side is v*t = t.
    # A point on a side can be defined by its distance from one vertex. Let's take the
    # counter-clockwise direction (V6->V1, V2->V3, V4->V5).
    # The initial distance from V6, V2, V4 is L/2 = 5.
    # The distance at time t is d(t) = 5 + t.
    # The parameter for linear interpolation P = A.lerp(B, lambda) is lambda = d(t)/L.
    lambda_t = (5 + t) / R

    # P(t) on side V6V1
    P_t = V6.lerp(V1, lambda_t)
    
    # Q(t) on side V2V3
    Q_t = V2.lerp(V3, lambda_t)
    
    # R(t) on side V4V5
    R_t = V4.lerp(V5, lambda_t)

    # Step 3: Calculate the square of the side length of the triangle T(t).
    # We calculate the squared distance between P(t) and Q(t).
    side_squared = P_t.distance(Q_t)**2
    
    # Sympy simplifies this expression.
    side_squared_simplified = sympy.simplify(side_squared)
    
    # Step 4: Calculate the area of the equilateral triangle.
    # The area is (sqrt(3)/4) * side^2.
    area_t = (sqrt3 / 4) * side_squared_simplified

    # Extract coefficients to print the formula clearly.
    poly_in_t = sympy.poly(side_squared_simplified, t)
    coeff_t2 = poly_in_t.coeff_monomial(t**2)
    constant_term = poly_in_t.coeff_monomial(1)

    # Step 5: Print the final result.
    print("The area of the triangle T(t) is A(t).")
    print("The squared side length s^2 of the triangle is a function of time t:")
    print(f"s^2(t) = {side_squared_simplified}")
    print("\nSince the triangle remains equilateral, its area A(t) = (sqrt(3)/4) * s^2(t).")
    print("The final expression for the area is:")
    
    s3_str = "sqrt(3)"
    denominator = 4
    
    print(f"\nA(t) = ({s3_str} / {denominator}) * ({coeff_t2}*t^2 + {constant_term})")

# Execute the calculation and print the result.
calculate_triangle_area()