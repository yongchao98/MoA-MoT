import sympy

def solve_particle_integral():
    """
    Solves the given physics problem by performing symbolic calculations.
    
    The steps are:
    1. Define the system matrix and its properties.
    2. Construct the general solution for the intermediate vector L(t).
    3. Use boundary conditions at t=tau to solve for integration constants.
    4. Determine the initial position vector r(0) as a function of tau.
    5. Calculate the required time-averaged integral.
    """
    
    # Define symbols
    t, tau = sympy.symbols('t tau', real=True, positive=True)
    c1, c2, c3 = sympy.symbols('c1 c2 c3')
    
    # Step 1: Define the system based on eigenvalues (0, 0, -4) and eigenvectors.
    # Eigenvector for lambda = 0
    v1 = sympy.Matrix([1, 0, -1])
    # Generalized eigenvector for lambda = 0
    u1 = sympy.Matrix([0, 1, -2])
    # Eigenvector for lambda = -4
    v3 = sympy.Matrix([3, 4, 5])
    
    # Step 2: Construct the general solution for L(t) = [Lx, Ly, Lz]
    # L(t) = (c1 + c2*t^2/2) * v1 + c2 * u1 + c3 * exp(-2*t^2) * v3
    L_t = (c1 + c2*t**2/2) * v1 + c2 * u1 + c3 * sympy.exp(-2*t**2) * v3
    
    # Step 3: Apply boundary conditions at t = tau.
    # At t=tau, r'(tau)=0, r(tau)=[0,0,1], which means L(tau) = [0, 0, -1].
    L_tau = L_t.subs(t, tau)
    
    # Set up the system of equations L(tau) = [0, 0, -1]
    # To make it easier for the solver, we can represent exp(-2*tau**2) as a variable
    # Let's solve it as a linear system for the terms.
    # eq1: c1 + c2*tau^2/2 + 3*c3*exp(-2*tau^2) = 0
    # eq2: c2 + 4*c3*exp(-2*tau^2) = 0
    # eq3: -(c1+c2*tau^2/2) - 2*c2 + 5*c3*exp(-2*tau^2) = -1
    # From eq2, c2 = -4*c3*exp(-2*tau^2).
    # Substitute into eq1, c1 = -c2*tau^2/2 - 3*c3*exp(-2*tau^2) = 2*c3*exp(-2*tau^2)*tau^2 - 3*c3*exp(-2*tau^2)
    # Substitute into eq3, we find 16*c3*exp(-2*tau^2) = -1.
    # So, c3*exp(-2*tau^2) = -1/16.
    
    sol = {
        c2: sympy.Rational(1, 4),
        c3: -sympy.exp(2*tau**2) / 16,
        c1: sympy.Rational(3, 16) - tau**2 / 8
    }

    # Substitute the solved constants back into L(t)
    L_t_solved = L_t.subs(sol)

    # Step 4: Find the initial position r(0) = -L(0).
    L_0 = L_t_solved.subs(t, 0)
    r_0 = -L_0

    # Calculate the sum S(tau) = x(0) + y(0) + z(0)
    # The sum is the sum of the elements of the r_0 vector.
    S_tau = r_0[0] + r_0[1] + r_0[2]
    S_tau_simplified = sympy.simplify(S_tau)

    # Step 5: Calculate the final integral
    integrand = 1 / S_tau_simplified
    integral_result = sympy.integrate(integrand, (tau, 0, sympy.oo))
    
    # Print the results
    print(f"The sum of initial positions S(τ) = x(0;τ) + y(0;τ) + z(0;τ) is found to be:")
    print(S_tau_simplified)
    print("\nThe integral to be calculated is ∫[0, ∞] 1/S(τ) dτ.")
    
    # To show the numbers in the final equation, we can express ln(16/9) as ln(16) - ln(9)
    # Extract numerator and denominator from the argument of the logarithm
    log_arg = sympy.exp(integral_result)
    num, den = sympy.fraction(log_arg)

    print(f"The result of the integration is ln({num}/{den}).")
    print(f"This can be written as ln({num}) - ln({den}).")
    print(f"The final numerical value is approximately: {integral_result.evalf()}")

solve_particle_integral()