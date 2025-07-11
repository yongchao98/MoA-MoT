import sympy
from sympy import sqrt, Rational

def solve_triangle_area():
    """
    This function calculates the area of the triangle T(t) as a function of time t.
    """
    # 1. Define symbols and constants from the problem
    t = sympy.Symbol('t')
    R = 10
    s = R  # Side length of the regular hexagon is equal to the radius
    v = 1  # Speed of the triangle's vertices

    # 2. Set up the hexagon's vertices in a fixed coordinate system
    # Using exact rational numbers for precision with sympy's symbolic computation
    A = sympy.Point(R, 0)
    B = sympy.Point(R * Rational(1, 2), R * sqrt(3) / Rational(2))
    C = sympy.Point(-R * Rational(1, 2), R * sqrt(3) / Rational(2))
    D = sympy.Point(-R, 0)
    E = sympy.Point(-R * Rational(1, 2), -R * sqrt(3) / Rational(2))
    F = sympy.Point(R * Rational(1, 2), -R * sqrt(3) / Rational(2))

    # 3. Parameterize the positions of the triangle's vertices P, Q, R
    # At t=0, the vertices are at the midpoints of sides AB, CD, EF.
    # A midpoint is at a distance of s/2 from a side's starting vertex.
    # The vertices move with speed v=1. We assume they move away from the
    # starting vertices (A, C, E) along their respective sides (AB, CD, EF).
    # The distance from the starting vertex at time t is:
    dist_from_start = s / Rational(2) + v * t

    # Calculate the coordinates of the triangle's vertices P(t), Q(t), R(t)
    # P is on side AB, starting from A
    P = A + (dist_from_start / s) * (B - A)
    # Q is on side CD, starting from C
    Q = C + (dist_from_start / s) * (D - C)
    # R is on side EF, starting from E
    R = E + (dist_from_start / s) * (F - E)

    # 4. Calculate the area using sympy's Polygon.area method (which uses the Shoelace formula)
    triangle = sympy.Polygon(P, Q, R)
    area_expr = sympy.simplify(triangle.area)

    # 5. Extract the numerical components of the formula to display them clearly.
    # The expected formula structure is: (factor * sqrt(root)) / denominator * (constant + t**2)
    
    # We isolate the rational part of the expression by dividing by sqrt(3)
    rational_part = sympy.expand(area_expr / sqrt(3))
    
    # Convert to a polynomial in t to safely extract coefficients
    poly_t = sympy.Poly(rational_part, t)
    
    # Get coefficients for t^2 and the constant term
    coeff_t2 = poly_t.coeff_monomial(t**2)
    coeff_t0 = poly_t.coeff_monomial(t**0)
    
    # From these coefficients, we derive the numbers for our desired formula structure
    factor = coeff_t2.p          # Numerator of the pre-factor (3)
    root = 3                     # The number inside the square root we divided out
    denominator = coeff_t2.q     # Denominator of the pre-factor (4)
    constant = coeff_t0 / coeff_t2 # The constant term inside the parenthesis (75)

    # Print the explanation and the final formula
    print("The area of a figure is not changed by rotation, so the angular velocity of the hexagon is irrelevant.")
    print("The area of the triangle T(t) is calculated based on the positions of its vertices as they move along the hexagon's sides.")
    print("\nThe derived formula for the area A(t) is:")
    print(f"A(t) = ({factor} * sqrt({root}) / {denominator}) * ({constant} + t^2)")
    print("\nEach number in the final equation is:")
    print(f"Factor: {factor}")
    print(f"Value inside square root: {root}")
    print(f"Denominator: {denominator}")
    print(f"Constant term in parenthesis: {constant}")

if __name__ == '__main__':
    solve_triangle_area()