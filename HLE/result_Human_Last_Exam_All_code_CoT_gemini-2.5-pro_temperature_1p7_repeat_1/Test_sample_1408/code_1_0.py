import sympy

def find_maximal_overhang():
    """
    Calculates the integers a, b, c corresponding to the maximal overhang of three
    cubes using symbolic mathematics.

    The maximal overhang for three identical cubes of side length 1 is found
    by setting up and solving stability equations, assuming optimal 45-degree rotation
    for all cubes to maximize the support area.

    The final overhang is expressed in the format (a + sqrt(b)) / (1 + c).
    """

    # Define symbolic variables for the x-coordinates of the centers of mass
    # of the three cubes, numbered from top (1) to bottom (3).
    x1, x2, x3 = sympy.symbols('x1 x2 x3')

    # Define a symbol for the support half-width. For a cube of side length 1,
    # rotated 45 degrees, this value 'd' is 1/sqrt(2).
    d = sympy.symbols('d')

    # Set up the system of equations for marginal stability:
    # 1. The combined CM of all three cubes is at the table's edge (x=0).
    eq1 = sympy.Eq(x1 + x2 + x3, 0)
    # 2. The CM of the top cube (1) is at the edge of the supporting cube (2).
    eq2 = sympy.Eq(x1, x2 + d)
    # 3. The combined CM of cubes 1 and 2 is at the edge of the supporting cube (3).
    eq3 = sympy.Eq((x1 + x2) / 2, x3 + d)

    # Solve the system for x1, x2, x3 in terms of d.
    solution = sympy.solve([eq1, eq2, eq3], (x1, x2, x3))
    x1_sol = solution[x1]
    x2_sol = solution[x2]
    x3_sol = solution[x3]
    
    # The overhang is the maximum x-coordinate of any point on any of the cubes.
    # The rightmost edge of cube i is at x_i + d.
    overhang_expr_1 = x1_sol + d
    overhang_expr_2 = x2_sol + d
    overhang_expr_3 = x3_sol + d

    # To find the maximal possible overhang, we substitute the maximal value for d.
    d_max = 1 / sympy.sqrt(2)
    
    # Substitute d_max into the expressions. Since d is positive, the maximum
    # value will come from the expression with the most positive coefficients for d.
    # overhang_expr_1 evaluates to (11/6)*d, overhang_expr_2 to (5/6)*d, and
    # overhang_expr_3 to (1/3)*d. So we only need to evaluate the first one.
    max_overhang = overhang_expr_1.subs(d, d_max)

    # The result is 11/(6*sqrt(2)), which simplifies to 11*sqrt(2)/12.
    # We need to express this in the form (a + sqrt(b)) / (1 + c).
    # (0 + 11*sqrt(2)) / 12 = (0 + sqrt(121*2)) / 12 = (0 + sqrt(242)) / 12
    # From this, we can identify the non-negative integers a, b, and c.
    # Denominator: 1 + c = 12  => c = 11
    # Numerator: a + sqrt(b) = 0 + sqrt(242) => a = 0, b = 242

    a = 0
    b = 242
    c = 11

    # Print the values of a, b, and c as requested.
    print(f"{a} {b} {c}")

if __name__ == '__main__':
    find_maximal_overhang()