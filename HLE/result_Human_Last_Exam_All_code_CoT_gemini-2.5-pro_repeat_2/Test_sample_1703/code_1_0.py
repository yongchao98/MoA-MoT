import sympy

def solve_triangle_area():
    """
    This function calculates the area of the triangle T(t) as a function of time t.
    """
    # 1. Define symbolic variable for time and constants
    t = sympy.Symbol('t')
    R = 10
    sqrt3 = sympy.sqrt(3)

    # 2. Define the hexagon's vertices in a stationary coordinate system.
    # The hexagon is centered at the origin, with one vertex on the positive x-axis.
    # Vertices are labeled A, B, C, D, E, F counter-clockwise.
    A = (R, 0)
    B = (R / 2, R * sqrt3 / 2)
    C = (-R / 2, R * sqrt3 / 2)
    D = (-R, 0)
    E = (-R / 2, -R * sqrt3 / 2)
    F = (R / 2, -R * sqrt3 / 2)

    # 3. Parameterize the positions of the triangle's vertices.
    # The triangle's vertices are on alternating sides, e.g., AB, CD, EF.
    # At t=0, they are at the midpoints. They move with speed v=1.
    # The side length of the hexagon is R.
    # The distance traveled from the start of the side is d(t) = R/2 + 1*t.
    # We use a parameter 'p' for interpolation along the side.
    p = (R / 2 + t) / R

    # P1 moves from A to B on side AB
    P1_x = (1 - p) * A[0] + p * B[0]
    P1_y = (1 - p) * A[1] + p * B[1]

    # P2 moves from C to D on side CD
    P2_x = (1 - p) * C[0] + p * D[0]
    P2_y = (1 - p) * C[1] + p * D[1]

    # P3 moves from E to F on side EF
    P3_x = (1 - p) * E[0] + p * F[0]
    P3_y = (1 - p) * E[1] + p * F[1]

    # 4. Calculate the area using the Shoelace formula.
    # Area = 0.5 * |(x1*y2 + x2*y3 + x3*y1) - (y1*x2 + y2*x3 + y3*x1)|
    area_expr = 0.5 * sympy.Abs(
        (P1_x * P2_y + P2_x * P3_y + P3_x * P1_y) -
        (P1_y * P2_x + P2_y * P3_x + P3_y * P1_x)
    )

    # 5. Simplify the area expression.
    simplified_area = sympy.simplify(area_expr)

    # 6. Extract the numbers from the final equation to display them clearly.
    # The expected form is (C2*t^2 + C0) * sqrt(3) / Denominator
    num, den = simplified_area.as_numer_denom()
    
    # The numerator should contain a sqrt(3) factor
    poly_in_t = sympy.simplify(num / sqrt3)

    # Extract coefficients from the polynomial part
    poly_in_t = sympy.poly(poly_in_t, t)
    c2 = poly_in_t.coeff_monomial(t**2)
    c0 = poly_in_t.coeff_monomial(1)

    # 7. Print the final result.
    print("The area of the triangle T(t) is independent of the hexagon's rotation.")
    print("The area A(t) as a function of time t is given by the formula:")
    print(f"A(t) = ({c2} * t^2 + {c0}) * sqrt(3) / {den}")
    print("\nWhere:")
    print(f"The coefficient of t^2 is {c2}.")
    print(f"The constant term is {c0}.")
    print(f"The expression is divided by {den}.")
    print("sqrt(3) is the square root of 3.")

# Execute the function to find and print the area formula.
solve_triangle_area()