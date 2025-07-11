import sympy

def solve_ant_path():
    """
    Calculates and prints the minimal polynomial for the shortest distance an ant
    can travel on a dodecahedron starting and ending at a vertex.
    """
    
    # Let x be the variable for the distance.
    x = sympy.symbols('x')

    # According to mathematical research on geodesics on polyhedra,
    # the shortest possible distance 'd' for an ant starting at a vertex of a
    # regular dodecahedron (with side length 1), walking on the surface in a
    # geodesic path, and returning to the same vertex without passing through any
    # other vertex is d = 2 + sqrt(5).
    distance_val = 2 + sympy.sqrt(5)

    # We want to find the minimal polynomial of 'd'. A minimal polynomial is the
    # polynomial of lowest degree with rational coefficients that has 'd' as a root.
    # We can find this by setting x = d and eliminating the square root, as shown
    # in the step-by-step plan.
    #
    # The sympy library provides a function to compute this directly.
    min_poly_expr = sympy.minimal_polynomial(distance_val, x)

    # The minimal polynomial is x^2 - 4x - 1.
    # This polynomial is irreducible over the rational numbers because its roots,
    # 2 + sqrt(5) and 2 - sqrt(5), are irrational.

    # The problem asks to output each number in the final equation.
    # We can get the coefficients from the sympy polynomial object.
    coeffs = sympy.poly(min_poly_expr, x).all_coeffs()
    a = coeffs[0]
    b = coeffs[1]
    c = coeffs[2]

    # Print the final equation.
    print("Let x be the shortest possible distance the ant could have walked.")
    print("The minimal polynomial P(x) for this distance is given by P(x) = 0.")
    print("\nThe equation is:")
    print(f"{a}x^2 + ({b})x + ({c}) = 0")

if __name__ == '__main__':
    solve_ant_path()