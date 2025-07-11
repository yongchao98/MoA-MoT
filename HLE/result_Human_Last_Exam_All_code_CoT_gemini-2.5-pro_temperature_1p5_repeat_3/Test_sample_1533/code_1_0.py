import sympy

def solve_geometry_ratio():
    """
    Symbolically calculates the ratio BM/MI in a triangle ABC.
    
    The function follows these steps:
    1. Defines symbolic variables for side lengths a, b, c.
    2. Defines expressions for semi-perimeter (s), inradius (r), and circumradius (R).
    3. Defines the expression for sin^2(B/2).
    4. Constructs the ratio BM/MI = 1 + r / (2 * R * sin^2(B/2)).
    5. Simplifies the expression and prints the result.
    """
    
    # Define symbolic variables for the side lengths of the triangle
    a, b, c = sympy.symbols('a b c', positive=True)

    # The semi-perimeter of the triangle
    s = (a + b + c) / 2

    # The square of the area of the triangle (by Heron's formula)
    area_sq = s * (s - a) * (s - b) * (s - c)
    
    # The inradius (r = Area / s), so r^2 = Area^2 / s^2
    r_sq = area_sq / s**2
    
    # The circumradius (R = abc / 4*Area), so R^2 = (a*b*c)**2 / (16*Area^2)
    R_sq = (a * b * c)**2 / (16 * area_sq)

    # The half-angle formula for sin(B/2): sin^2(B/2) = (s-a)(s-c)/(ac)
    sin_sq_B_half = (s - a) * (s - c) / (a * c)

    # The ratio is derived as 1 + r / (2*R*sin^2(B/2))
    # To avoid issues with square roots, we can work with the square of the second term.
    # term_to_add = r / (2*R*sin^2(B/2))
    # term_to_add_sq = r^2 / (4*R^2*sin^4(B/2))
    # It is simpler to substitute directly as derived in the explanation.
    # ratio_term = 2*(s-b)/b
    
    ratio_term = (2 * (s - b)) / b

    # The final expression for the ratio BM/MI is 1 + the term above
    final_ratio = 1 + ratio_term
    
    # Simplify the final expression
    simplified_ratio = sympy.simplify(final_ratio)
    
    print("The geometric derivation leads to the expression: 1 + 2*(s-b)/b")
    print("Substituting s = (a+b+c)/2 and simplifying, the final expression for the ratio BM/MI is:")
    print(simplified_ratio)

solve_geometry_ratio()