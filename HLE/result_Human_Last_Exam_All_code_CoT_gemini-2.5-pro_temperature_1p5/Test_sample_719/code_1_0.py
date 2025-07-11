import sympy as sp

def symbolic_derivation():
    """
    This function performs the symbolic derivation to find theta'(t) as described in the problem.
    """
    # Step 1: Define all symbolic variables and functions
    t = sp.Symbol('t')
    f = sp.Function('f')(t)
    theta = sp.Function('theta')(t)
    r = sp.Function('r')(t)

    # Define their derivatives
    f_prime = f.diff(t)
    r_prime = r.diff(t)
    theta_prime = theta.diff(t)

    # Step 2: Set up the system of equations.
    # This system results from equating two different ways of expressing the time-derivatives of the components (w1, w2).
    # First way is from the ODE system derived for w1, w2 when K=0:
    #   w1_dot = -(f'/f) * w1
    #   w2_dot = f * w1
    # Second way is by differentiating the polar coordinate representation w1=r*cos(theta), w2=r*sin(theta).
    # This leads to the following system of linear equations for r_prime and theta_prime.
    # We can divide by r (assuming a non-trivial solution, r != 0).
    eq1 = sp.Eq(r_prime * sp.cos(theta) - r * theta_prime * sp.sin(theta), -f_prime/f * r * sp.cos(theta))
    eq2 = sp.Eq(r_prime * sp.sin(theta) + r * theta_prime * sp.cos(theta), f * r * sp.cos(theta))

    # Step 3: Solve the system for theta_prime
    # The `solve` function in SymPy can handle this system of two equations and two unknowns (r_prime, theta_prime).
    solution = sp.solve([eq1, eq2], (r_prime, theta_prime), dict=True)
    
    # Extract the solution for theta_prime.
    theta_prime_solution = solution[0][theta_prime]

    # Step 4: Print the resulting expression.
    # The result has two main terms.
    # Term 1: f(t) * cos(theta(t))**2
    # Term 2: (f'(t)/f(t)) * sin(theta(t)) * cos(theta(t))
    # We print the full simplified expression for theta'(t).
    final_expression = sp.simplify(theta_prime_solution)

    print("The final equation for theta'(t) is:")
    print(final_expression)

# Execute the symbolic derivation
symbolic_derivation()