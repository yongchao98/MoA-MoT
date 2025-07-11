import sympy

def solve_srpt_time():
    """
    Calculates the optimal mean response time x for an M/G/1 queue
    with lambda=3/2 and service time U(0,1) under SRPT policy.
    Then, it identifies the remaining term after removals.
    """
    # Define symbols and constants
    s = sympy.Symbol('s')
    lambda_rate = sympy.Rational(3, 2)
    
    # PDF f(s) is 1 on [0, 1]
    
    # Calculate rho_s and E[S_s^2]
    # In general, rho_s = lambda_rate * integrate(y*f(y), (y, 0, s))
    # E[S_s^2] = integrate(y**2*f(y), (y, 0, s))
    # For f(s)=1 on [0,1], these are:
    rho_s = lambda_rate * (s**2 / 2)
    E_Ss_sq = s**3 / 3
    
    # Define the two parts of the E[W] integrand
    integrand_part1 = (lambda_rate / 2) * E_Ss_sq / (1 - rho_s)**2
    integrand_part2 = s * rho_s / (1 - rho_s)
    
    # Calculate the two integrals for E[W]
    integral1 = sympy.integrate(integrand_part1, (s, 0, 1))
    integral2 = sympy.integrate(integrand_part2, (s, 0, 1))
    
    # Calculate E[W]
    E_W = integral1 + integral2
    
    # Mean service time E[S] for U(0,1) is 1/2
    E_S = sympy.Rational(1, 2)
    
    # Calculate optimal mean response time x = E[T]
    x = E_S + E_W
    
    # The result is x = 2/3 + (8/9)*ln(2)
    # Extract the coefficients
    # The rational part is the constant term in the expression for x.
    # The Add constructor separates the terms.
    rational_term = sympy.S(0)
    log_term = sympy.S(0)

    # Decompose the expression x
    if x.is_Add:
        for term in x.args:
            if term.is_rational:
                rational_term += term
            else:
                log_term += term
    elif x.is_rational:
        rational_term = x
    else:
        log_term = x

    # The log term is (8/9)*ln(2). We extract its coefficient and argument.
    log_coeff = log_term.as_coeff_Mul()[0]
    log_base = log_term.as_coeff_Mul()[1].args[0]
    
    print("The optimal mean response time x is given by the equation:")
    print(f"x = {rational_term} + {log_coeff}*ln({log_base})")
    
    print("\nTask: Remove all additive rational terms and all additive terms which are logarithms of rational numbers from x.")
    print(f"1. The term '{rational_term}' is a rational term. It is removed.")
    print(f"2. The term '{log_term}' is not rational. It can be written as ln({log_base}**{log_coeff}), which is not a logarithm of a rational number. It remains.")
    
    remaining_term_coeff = log_coeff
    remaining_term_log_arg = log_base
    
    # Formatting for LaTeX
    # Numerator and denominator of the rational coefficient
    num = remaining_term_coeff.p
    den = remaining_term_coeff.q
    latex_string = f"\\frac{{{num}}}{{{den}}}\\ln{{{remaining_term_log_arg}}}"
    
    print(f"\nThe remaining term is {log_coeff}*ln({log_base}).")
    print(f"Formatted in LaTeX: {latex_string}")

solve_srpt_time()