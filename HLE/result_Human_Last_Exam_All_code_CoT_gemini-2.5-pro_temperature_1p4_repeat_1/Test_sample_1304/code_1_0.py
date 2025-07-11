import sympy

def solve_max_c3():
    """
    This function solves the maximization problem to find the maximum value of c_3.

    The problem reduces to maximizing the function c_3 = -35/2 * z1*z2*(z1+z2)
    subject to -1 <= z1 < 0 and 0 < z2 <= 1.
    By substituting x = -z1 and y = z2, the problem becomes maximizing
    G(x, y) = 35/2 * x*y*(y-x) for x, y in (0, 1].

    The analytical derivation shows that the maximum occurs on the boundary where y=1.
    We use sympy to find this maximum.
    """

    # Define symbolic variables
    x, y = sympy.symbols('x y')

    # Define the function G(x, y) which is proportional to c_3
    # G = x * y * (y - x) represents the variable part of the expression for c_3
    G = x * y * (y - x)

    # We analyze the boundary where y = 1 and x is in (0, 1].
    G_y1 = G.subs(y, 1)

    # To find the maximum of G(x, 1), we take its derivative with respect to x.
    dG_dx = sympy.diff(G_y1, x)

    # Solve for x where the derivative is zero to find the critical point.
    critical_points = sympy.solve(dG_dx, x)
    x_max_val = critical_points[0]

    # Substitute this value of x back into G(x, 1) to find its maximum value.
    max_G_val = G_y1.subs(x, x_max_val)

    # The maximum value of c_3 is (35/2) * max_G_val.
    c3_factor_num = 35
    c3_factor_den = 2

    final_c3_val = (sympy.S(c3_factor_num)/c3_factor_den) * max_G_val

    # Print the steps of the final calculation
    max_G_num = sympy.fraction(max_G_val)[0]
    max_G_den = sympy.fraction(max_G_val)[1]
    final_c3_num = sympy.fraction(final_c3_val)[0]
    final_c3_den = sympy.fraction(final_c3_val)[1]
    
    print(f"The maximization problem reduces to finding the maximum of a function of two variables.")
    print(f"Analysis shows the maximum is on a boundary, leading to maximizing G(x) = x*(1-x).")
    print(f"The maximum value of G(x) is {max_G_val}, found at x = {x_max_val}.")
    print(f"The maximum value of c_3 is obtained from this result.")
    print(f"Final equation:")
    print(f"{c3_factor_num} / {c3_factor_den} * {max_G_num} / {max_G_den} = {final_c3_num} / {final_c3_den}")

if __name__ == '__main__':
    solve_max_c3()