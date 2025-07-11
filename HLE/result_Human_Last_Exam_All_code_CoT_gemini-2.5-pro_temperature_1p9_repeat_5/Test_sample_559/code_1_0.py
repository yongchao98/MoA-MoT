import sympy

def find_separatrix_equation():
    """
    Finds the separatrix for the given system of ODEs by assuming a polynomial form
    d(u) = a*u^2 + b*u + c, and solving for its coefficients using symbolic math.
    """
    # Define symbolic variables
    u = sympy.Symbol('u')
    a, b, c = sympy.symbols('a b c')

    # Define the hypothesized polynomial solution and its derivative
    d_poly = a * u**2 + b * u + c
    d_prime = sympy.diff(d_poly, u)

    # The original ODE is: d'/u' = N/D, which we rewrite as D*d' - N = 0 to avoid fractions.
    # d'(t) = 2*d**2 + (-3*u+5*u**2)*d - u**3*(1-u)
    # u'(t) = u**2*(u-1) = u**3 - u**2
    # So the phase plane ODE is: (u**3 - u**2)*d' = 2*d**2 + (5*u**2 - 3*u)*d + u**4 - u**3
    
    # Left hand side of the equation
    lhs = (u**3 - u**2) * d_prime
    # Right hand side of the equation
    rhs = 2 * d_poly**2 + (5 * u**2 - 3 * u) * d_poly + u**4 - u**3

    # We need to solve lhs - rhs = 0.
    # Substitute the polynomial and simplify the expression.
    equation = sympy.simplify(lhs - rhs)

    # The 'equation' is a polynomial in u. For it to be zero for all u,
    # all its coefficients must be zero. We extract these coefficients.
    poly_equation = sympy.Poly(equation, u)
    coeffs = poly_equation.all_coeffs()

    # Solve the system of algebraic equations for the coefficients (a, b, c).
    solutions = sympy.solve(coeffs, (a, b, c), dict=True)

    # From our analysis, the saddle point is at (d, u) = (-1, 1).
    # The separatrix must pass through this point. We test our solutions.
    saddle_u = 1
    saddle_d = -1
    
    print("Found potential quadratic solutions d(u) = a*u^2 + b*u + c. The coefficients (a,b,c) are:")
    print(solutions)
    
    final_a, final_b, final_c = None, None, None

    for sol in solutions:
        a_sol, b_sol, c_sol = sol.get(a), sol.get(b), sol.get(c)
        # Construct the polynomial for this solution
        d_sol_expr = a_sol * u**2 + b_sol * u + c_sol
        
        # Check if it passes through the saddle point d(1) = -1
        if sympy.Eq(d_sol_expr.subs(u, saddle_u), saddle_d):
            print(f"\nChecking solution d(u) = {d_sol_expr}...")
            print(f"This solution passes through the saddle point ({saddle_d}, {saddle_u}). This is the separatrix.")
            final_a = a_sol
            final_b = b_sol
            final_c = c_sol
            break
            
    if final_a is not None:
        print("\nThe equation for the separatrix is d = a*u^2 + b*u + c")
        print("with the following coefficients:")
        print(f"a = {final_a}")
        print(f"b = {final_b}")
        print(f"c = {final_c}")
        
        print("\nFinal equation is:")
        # Print each number in the final equation.
        print(f"d = ({final_a}) * u^2 + ({final_b}) * u + ({final_c})")
    else:
        print("No quadratic polynomial solution found that acts as the separatrix.")

find_separatrix_equation()