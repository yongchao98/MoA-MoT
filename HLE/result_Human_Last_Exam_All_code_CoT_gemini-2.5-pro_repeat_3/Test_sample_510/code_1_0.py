import sympy

def solve_queueing_problem():
    """
    This function calculates the optimal mean response time 'x' for the given
    queueing system and isolates the term specified in the problem.
    """
    # Step 1 & 2: Define symbols and parameters
    s, t = sympy.symbols('s t')
    lambda_rate = sympy.Rational(3, 2)
    
    # For U(0, 1), f(s) = 1 and F(s) = s on the interval [0, 1]
    f_s = 1
    F_s = s
    
    # Calculate rho(s) = lambda * integral from 0 to s of (t * f(t)) dt
    rho_s = lambda_rate * sympy.integrate(t * f_s, (t, 0, s))
    
    # The integrand is (1 - F(s)) / (1 - rho(s))
    integrand = (1 - F_s) / (1 - rho_s)
    
    # Step 3: Calculate the integral for x from 0 to 1
    # For s > 1, F(s) = 1, so the integrand is 0.
    x = sympy.integrate(integrand, (s, 0, 1))
    
    # The integral can be manually split for clarity:
    # integrand = (1-s)/(1-3*s**2/4) = (4-4s)/(4-3s**2)
    # This splits into integral of 4/(4-3s**2) and -4s/(4-3s**2)
    part1_integrand = 4 / (4 - 3*s**2)
    part2_integrand = -4*s / (4 - 3*s**2)
    
    term1 = sympy.integrate(part1_integrand, (s, 0, 1))
    term2 = sympy.integrate(part2_integrand, (s, 0, 1))

    # Step 4: Analyze the result
    # x is the sum of term1 and term2.
    # term1 is the irrational/transcendental part.
    # term2 is the logarithm of a rational number part.
    
    print("The optimal mean response time x is given by the integral:")
    print(f"x = Integral from 0 to 1 of ({integrand}) ds")
    print("\nThe integral evaluates to:")
    print(f"x = {sympy.pretty(x)}")
    print("\nThis can be broken down into two terms:")
    print(f"Term A = {sympy.pretty(term1)}")
    print(f"Term B = {sympy.pretty(term2)}")

    print("\nAccording to the problem, we must remove:")
    print("1. Additive rational terms (there are none).")
    print("2. Additive terms which are logarithms of rational numbers.")
    print(f"Term B, which is {sympy.pretty(term2)}, is a rational multiple of a logarithm of a rational number (2), so it is removed.")
    
    # Step 5: Identify the remaining term
    remaining_term = term1
    print("\nThe remaining term of x is:")
    print(sympy.pretty(remaining_term))
    
    # Format for LaTeX output
    # Rational: 2/3, Algebraic irrational: sqrt(3), Transcendental: ln(2+sqrt(3))
    latex_answer = sympy.latex(sympy.Rational(2,3)) + r' \sqrt{3} \ln(2+\sqrt{3})'
    print(f"\nFinal Answer in LaTeX format: {latex_answer}")


solve_queueing_problem()
