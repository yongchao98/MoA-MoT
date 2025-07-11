import sympy

def solve_and_format():
    """
    This function uses sympy to symbolically solve for the optimal mean response time x,
    and then isolates the remaining term based on the problem's instructions.
    """
    # Define symbols and parameters
    s, y = sympy.symbols('s y')
    lambda_val = sympy.Rational(3, 2)

    # Calculate rho(s) for U(0,1) distribution (f(y)=1)
    rho_s = lambda_val * sympy.integrate(y, (y, 0, s))

    # Define the integrand for the mean response time x
    # E[T_s] = s/(1-rho(s)) + (lambda * integral(y^2 dy)) / (2*(1-rho(s))^2)
    numerator_term2 = lambda_val * sympy.integrate(y**2, (y, 0, s))
    denominator_term2 = 2 * (1 - rho_s)**2
    term1_integrand = s / (1 - rho_s)
    term2_integrand = numerator_term2 / denominator_term2
    integrand = term1_integrand + term2_integrand

    # Integrate E[T_s]*f(s) from 0 to 1. Since f(s)=1, just integrate.
    x = sympy.integrate(integrand, (s, 0, 1))

    # Now, analyze the expression for x
    # x is expected to be an Add object, e.g., 2/3 + 8*ln(2)/9
    if isinstance(x, sympy.Add):
        terms = x.args
    else:
        terms = [x]

    remaining_term = sympy.S.Zero

    for term in terms:
        # Check if term is rational
        if term.is_rational:
            continue # Remove rational terms

        # Check if term is a logarithm of a rational number, i.e., term = ln(q) for q in Q.
        # This is true if exp(term) is rational.
        # For a term c*ln(r) with r rational, this holds if r^c is rational.
        # For prime r, this requires c to be an integer.
        is_log_of_rational = False
        # Deconstruct term as coeff*log(arg)
        coeff, log_part = term.as_coeff_log()
        if log_part is not None:
            # Check if log argument is rational and coefficient is integer
            if log_part.args[0].is_rational and coeff.is_integer:
                is_log_of_rational = True
        
        if not is_log_of_rational:
            remaining_term += term
    
    # Format the remaining term into LaTeX
    # For 8/9*ln(2), coeff=8/9, terms=(ln(2),)
    coeff, non_rat_parts = remaining_term.as_coeff_mul(rational=True)
    
    # We will print the numbers in the final expression as requested.
    # For a/b * ln(c), the numbers are a, b, c.
    print(f"The optimal mean response time is x = {x}.")
    print(f"After removing rational terms and logarithms of rational numbers, the remaining term is: {remaining_term}")
    
    # Deconstruct for printing numbers
    a = coeff.p
    b = coeff.q
    # Assume the non-rational part is a single log term
    log_arg = non_rat_parts[0].args[0]
    c = log_arg
    
    print("The final expression for the remaining term is composed of the following numbers:")
    print(f"Rational multiplicand: {a}/{b}")
    print(f"Argument of the logarithm: {c}")
    
    final_latex_str = sympy.latex(remaining_term)
    print(f"\nFinal Answer in LaTeX format: {final_latex_str}")

solve_and_format()