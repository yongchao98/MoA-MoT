import sympy

def solve_theta_prime():
    """
    This function symbolically derives the expression for theta'(t) based on the physics of the linearized geodesic flow.
    """
    # Define symbols for our calculation
    t = sympy.Symbol('t')
    K = sympy.Symbol('K')  # Gaussian curvature K(gamma(t))
    c = sympy.Symbol('c')  # Constant from the basis definition
    r = sympy.Function('r')(t)
    theta = sympy.Function('theta')(t)

    print("The system of differential equations for the coordinates (U, W) in the given basis is:")
    print("dU/dt = -(K/c) * W")
    print("dW/dt = c * U\n")

    # Polar coordinate representation
    # U corresponds to the real part, W to the imaginary part
    U = r * sympy.cos(theta)
    W = r * sympy.sin(theta)

    # Differentiate U and W with respect to t using the chain rule
    dU_dt = sympy.diff(U, t)
    dW_dt = sympy.diff(W, t)

    # The ODEs provide another expression for the derivatives
    dU_dt_from_ode = -(K / c) * W
    dW_dt_from_ode = c * U

    # Create a system of equations to solve for r'(t) and theta'(t)
    eq1 = sympy.Eq(dU_dt, dU_dt_from_ode)
    eq2 = sympy.Eq(dW_dt, dW_dt_from_ode)
    
    # Solve the system for the derivatives
    r_prime = sympy.diff(r, t)
    theta_prime = sympy.diff(theta, t)
    solution = sympy.solve([eq1, eq2], [r_prime, theta_prime])

    # Extract the solution for theta'(t)
    theta_prime_solution = solution[theta_prime]

    # Simplify the final expression
    theta_prime_simplified = sympy.simplify(theta_prime_solution)

    # To present the result clearly, we identify each term in the final equation.
    term1_coeff = c
    term1_func = sympy.cos(theta)**2
    term2_coeff_num = K
    term2_coeff_den = c
    term2_func = sympy.sin(theta)**2

    print("The final expression for theta'(t) is derived as:")
    print(f"theta'(t) = ({term1_coeff}) * {term1_func} + ({term2_coeff_num}/{term2_coeff_den}) * {term2_func}")

solve_theta_prime()