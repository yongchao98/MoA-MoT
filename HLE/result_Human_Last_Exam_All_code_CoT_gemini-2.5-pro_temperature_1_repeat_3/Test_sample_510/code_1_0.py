import sympy

def solve_queueing_problem():
    """
    Calculates the optimal mean response time for the given M/G/1 system,
    isolates the specified term, and prints it in LaTeX format.
    """
    # Step 1 & 2: Define symbols and parameters
    s, y = sympy.symbols('s y')
    lamda = sympy.Rational(3, 2)
    f_s = 1  # PDF of U(0, 1) job size distribution

    # Step 3: Calculate the components for the SRPT mean response time formula
    
    # Mean job size, E[S]
    E_S = sympy.integrate(s * f_s, (s, 0, 1))
    
    # Overall traffic intensity, rho = lambda * E[S]
    rho = lamda * E_S
    
    # s-dependent traffic intensity, rho_s = lambda * integral from 0 to s of y*f(y)dy
    rho_s = lamda * sympy.integrate(y * f_s, (y, 0, s))
    
    # The mean waiting time E[W] is the sum of two integrals, I1 and I2.
    
    # First integral part of E[W]
    # Integrand for I1
    term1_numerator = lamda * sympy.integrate(y**2 * f_s, (y, 0, s))
    term1_denominator = 2 * (1 - rho_s)**2
    term1_integrand = (term1_numerator / term1_denominator) * f_s
    I1 = sympy.integrate(term1_integrand, (s, 0, 1))

    # Second integral part of E[W]
    # Integrand for I2
    term2_integrand = ((rho - rho_s) / (1 - rho_s)) * s * f_s
    I2 = sympy.integrate(term2_integrand, (s, 0, 1))
    
    # Total mean waiting time E[W]
    E_W = I1 + I2
    
    # Optimal mean response time, x = E[T] = E[S] + E[W]
    x = E_S + E_W

    # Step 4: Isolate the required term from x
    
    # Expand the expression for x to separate the additive terms
    x_expanded = sympy.expand(x)
    
    remaining_term = sympy.S(0)
    
    # The additive terms of x
    additive_terms = x_expanded.as_ordered_terms()
    
    for term in additive_terms:
        # Condition 1: Remove rational terms
        if term.is_rational:
            continue
        # Condition 2: Remove terms which are logarithms of rational numbers.
        # A term 't' is a log of a rational if exp(t) is rational.
        elif sympy.exp(term).is_rational:
            continue
        else:
            # This term is not removed
            remaining_term += term
            
    # Step 5: Format the result in LaTeX
    # Use ln_notation for natural log and clean up the output for readability
    latex_output = sympy.latex(remaining_term, ln_notation=True)
    latex_output = latex_output.replace(r'\left(', '').replace(r'\right)', '')
    latex_output = latex_output.replace('{', ' ').replace('}', ' ')

    print("The optimal mean response time x is:")
    sympy.pretty_print(x)
    print("\nAfter removing rational terms and logarithms of rational numbers, the remaining term is:")
    print(latex_output)

solve_queueing_problem()