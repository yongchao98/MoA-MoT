import sympy

def find_separatrix():
    """
    Finds the separatrix for the given system of differential equations.
    """
    # Define symbols
    d, u = sympy.symbols('d u')

    # Define the system of ODEs
    # d'(t) = 2*d**2 + (-3*u + 5*u**2)*d - u**3*(1-u)
    # u'(t) = (u-1)*u**2
    d_prime = 2*d**2 + (-3*u + 5*u**2)*d - u**2*u*(1-u)
    u_prime = (u-1)*u**2

    print("--- Step 1: Finding Equilibrium Points ---")
    u_sols = sympy.solve(u_prime, u)
    equilibria = []
    for u_sol in u_sols:
        # We are only interested in the u >= 0 half-plane
        if u_sol.is_real and u_sol >= 0:
            d_eq = d_prime.subs(u, u_sol)
            d_sols = sympy.solve(d_eq, d)
            for d_sol in d_sols:
                equilibria.append((d_sol, u_sol))

    print(f"The equilibrium points (d, u) in the u>=0 plane are: {equilibria}")
    
    print("\n--- Step 2: Stability Analysis via Jacobian ---")
    system_matrix = sympy.Matrix([d_prime, u_prime])
    vars_matrix = sympy.Matrix([d, u])
    J = system_matrix.jacobian(vars_matrix)

    saddle_points = []
    for point in equilibria:
        J_at_point = J.subs([(d, point[0]), (u, point[1])])
        eigenvals = list(J_at_point.eigenvals().keys())
        
        # Check if eigenvalues have opposite signs (saddle point condition)
        is_saddle = len(eigenvals) == 2 and eigenvals[0] * eigenvals[1] < 0
        point_type = "a saddle point" if is_saddle else "not a saddle point"
        
        print(f"Point {point} is {point_type} with eigenvalues {eigenvals}.")
        if is_saddle:
            saddle_points.append(point)

    print(f"\nThe saddle point is: {saddle_points[0]}")
    
    print("\n--- Step 3: Finding the Invariant Manifold ---")
    print("A separatrix is an invariant manifold, often related to a saddle point.")
    print("Let's test for a solution of the form d = c * u**2.")

    c = sympy.Symbol('c')
    # Proposed solution d = c * u^2
    # Its time derivative is d'(t) = d(c*u^2)/du * u'(t) = 2*c*u * u'(t)
    d_test = c*u**2
    d_prime_test = sympy.diff(d_test, u) * u_prime
    
    # Substitute d = c*u^2 into the original equation for d'(t)
    d_prime_orig = d_prime.subs(d, d_test)
    
    # The equation must hold for all u, so the difference must be zero.
    equation_for_c = sympy.simplify(d_prime_orig - d_prime_test)
    
    # Since this polynomial in u must be zero for all u, all its coefficients must be zero.
    # Sympy can solve for c directly.
    sol = sympy.solve(equation_for_c, c)
    c_val = sol[0]
    
    print(f"Substituting d = c*u^2 into the system yields the equation: {equation_for_c} = 0")
    print(f"Solving for the constant c gives c = {c_val}.")
    
    # Construct the final separatrix equation: d - c*u^2 = 0
    separatrix_eq = d - c_val*u**2
    
    print("\n--- Step 4: The Separatrix Equation ---")
    print(f"The separatrix is given by the equation: {separatrix_eq} = 0")
    
    print("\nThe equation for the separatrix with explicit coefficients is:")
    d_coeff = separatrix_eq.coeff(d)
    u2_coeff = separatrix_eq.coeff(u**2)
    # The term u**2 has exponent 2.
    print(f"{d_coeff}*d(t) + {u2_coeff}*u(t)^2 = 0")


if __name__ == '__main__':
    find_separatrix()