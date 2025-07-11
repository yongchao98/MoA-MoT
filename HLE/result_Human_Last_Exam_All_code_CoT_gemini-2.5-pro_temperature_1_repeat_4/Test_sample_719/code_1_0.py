import sympy as sp

def solve_for_theta_prime():
    """
    This function uses the SymPy library to symbolically derive the expression for theta'(t).
    """
    # Define time and the functions of time
    t = sp.Symbol('t')
    r = sp.Function('r')(t)
    theta = sp.Function('theta')(t)
    f = sp.Function('f')(t)
    y = sp.Function('y')(t)

    # Define the system of equations based on the physical model
    # Eq1: y(t) = r(t)*sin(theta(t))
    # Eq2: y_dot(t) = r(t)*cos(theta(t))*f(t)
    # Eq3: y_ddot(t) = 0
    eq1_rhs = r * sp.sin(theta)
    eq2_rhs = r * sp.cos(theta) * f

    # From Eq1 and Eq2, we can derive an expression for r_dot = dr/dt
    # Differentiate Eq1 with respect to t
    y_dot_from_eq1 = sp.diff(eq1_rhs, t)
    
    # Equate this with the RHS of Eq2 to form an equation for r_dot
    r_dot_eq = sp.Eq(y_dot_from_eq1, eq2_rhs)

    # Solve for r_dot
    r_dot = sp.diff(r, t)
    r_dot_sol = sp.solve(r_dot_eq, r_dot)[0]

    # Now, use the condition y_ddot = 0. Differentiate the RHS of Eq2.
    y_ddot_expr = sp.diff(eq2_rhs, t)

    # Substitute the solution for r_dot into the y_ddot expression
    y_ddot_final_expr = y_ddot_expr.subs(r_dot, r_dot_sol)

    # The result must be zero. We solve for theta_dot = d(theta)/dt.
    theta_dot = sp.diff(theta, t)
    theta_dot_eq = sp.Eq(y_ddot_final_expr, 0)
    
    # Solve for theta_dot
    theta_dot_sol = sp.solve(theta_dot_eq, theta_dot)[0]

    # Print the result
    print("The derived expression for theta'(t) is:")
    sp.pprint(theta_dot_sol)

    # Show the final equation with explicit numbers as requested
    f_prime = sp.Derivative(f, t)
    final_eq_str = f"theta'(t) = (1) * {f}*cos({theta})**2 + (1) * ({f_prime}/{f})*cos({theta})*sin({theta})"
    
    print("\nThis corresponds to option F.")
    print("\nThe final equation with coefficients displayed is:")
    print(final_eq_str)


if __name__ == '__main__':
    solve_for_theta_prime()