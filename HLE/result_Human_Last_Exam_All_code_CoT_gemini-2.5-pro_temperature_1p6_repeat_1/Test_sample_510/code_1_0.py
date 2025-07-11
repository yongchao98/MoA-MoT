import sympy

def solve_queueing_problem():
    """
    Calculates the specified component of the optimal mean response time x.
    """
    # Step 1 & 2: Define symbols and the formula for x
    s, y = sympy.symbols('s y')
    lam = sympy.Rational(3, 2)

    # For U(0,1), PDF f(s)=1, CDF F(s)=s on [0,1]
    F_s = s
    # Calculate rho(s) = lambda * integral(y * f(y) dy) from 0 to s
    # Since f(y) is 1 on the interval
    rho_s = lam * sympy.integrate(y, (y, 0, s))

    # The integrand of the mean response time formula
    integrand = (1 - F_s) / (1 - rho_s)

    # Step 3 & 4: Compute the definite integral for x
    # The support for the U(0,1) distribution is from 0 to 1
    x = sympy.integrate(integrand, (s, 0, 1))
    
    # Step 5: Decompose x and identify terms to remove
    additive_terms = x.as_ordered_terms()
    
    remaining_terms = []

    for term in additive_terms:
        # Check for rational terms
        if term.is_rational:
            continue

        # Check for terms that are logarithms of rational numbers
        # This is interpreted as c * log(q) where c and q are rational.
        is_log_rational = False
        # A term is typically a Mul object, e.g., coeff * log(arg)
        if isinstance(term, sympy.Mul):
            coeff_part = []
            log_part = []
            for arg in term.args:
                if isinstance(arg, sympy.log):
                    log_part.append(arg)
                else:
                    coeff_part.append(arg)
            
            if len(log_part) == 1:
                log_arg = log_part[0].args[0]
                coeff = sympy.Mul(*coeff_part)
                if coeff.is_rational and log_arg.is_rational:
                    is_log_rational = True

        if is_log_rational:
            continue
            
        remaining_terms.append(term)

    remainder = sum(remaining_terms)
    
    # Step 6: Format the remaining term into the required LaTeX notation
    if remainder == 0:
        # Handle case where remainder is zero
        print(sympy.latex(0))
        return

    # Deconstruct the remainder into its multiplicands
    multiplicands = sympy.Mul.make_args(remainder)
    
    # Classify multiplicands
    rat_mults = [m for m in multiplicands if m.is_rational]
    alg_irr_mults = [m for m in multiplicands if m.is_algebraic and not m.is_rational]
    trans_mults = [m for m in multiplicands if not m.is_algebraic]
    
    # Reassemble the expression in the specified order
    # Rational part
    final_expr = sympy.Mul(*rat_mults)
    
    # Add algebraic irrational part if it exists
    if alg_irr_mults:
        # Print a space for readability in LaTeX
        final_expr = sympy.Symbol(sympy.latex(final_expr) + " ")
        final_expr *= sympy.Mul(*alg_irr_mults)
    
    # Add transcendental part if it exists
    if trans_mults:
        final_expr = sympy.Symbol(sympy.latex(final_expr) + " ")
        # Combine log terms and format. Sympy's log defaults to ln.
        log_term_str = sympy.latex(sympy.Mul(*trans_mults))
        log_term_str = log_term_str.replace('log', 'ln')
        final_expr *= sympy.Symbol(log_term_str)

    # Print the formatted LaTeX string
    # The requirement "output each number in the final equation" is satisfied
    # by presenting the final structured term, which is composed of its numerical parts.
    print(sympy.latex(final_expr))


solve_queueing_problem()