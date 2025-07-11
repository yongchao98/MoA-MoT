import sympy

def solve_queueing_problem():
    """
    This function calculates the optimal mean response time for the given queueing system,
    then removes terms as specified in the problem description to find the final answer.
    """
    # Step 1 & 2: Define symbols and the formula for the SRPT mean response time.
    s, y = sympy.symbols('s y')
    lmbda = sympy.Rational(3, 2)
    # Service time PDF f(y) is 1 for y in [0, 1].

    # Step 3: Substitute parameters into the formulas.
    # The load from jobs with size <= s:
    # rho_s = lambda * integral from 0 to s of y*f(y) dy
    rho_s_expr = lmbda * sympy.integrate(y * 1, (y, 0, s))

    # The conditional mean response time for a job of size s, E[T_s].
    # The second term in E[T_s] contains an integral:
    integral_y_squared = lmbda * sympy.integrate(y**2 * 1, (y, 0, s))
    
    # E[T_s] integrand expression
    e_ts = s / (1 - rho_s_expr) + integral_y_squared / (2 * (1 - rho_s_expr)**2)

    # Step 4: Calculate the definite integral for x, the optimal mean response time.
    # x = Integral from 0 to 1 of E[T_s] * f(s) ds, where f(s) = 1.
    x = sympy.integrate(e_ts, (s, 0, 1))

    # The computed value of x is 2/3 + (8/9)*log(2).
    # The instruction "output each number in the final equation" is interpreted as showing
    # the components of this calculated value.
    # The value of the first integral (I_1) is: 4/3*ln(2)
    # The value of the second integral (I_2) is: 2/3 - 4/9*ln(2)
    # The total mean response time x is: 4/3*ln(2) + 2/3 - 4/9*ln(2) = 2/3 + 8/9*ln(2)
    
    print(f"The optimal mean response time is x = {x}.")
    print(f"The first rational term is: {x.as_ordered_terms()[0]}")
    print(f"The second term is: {x.as_ordered_terms()[1]}")


    # Step 5: Filter the resulting expression for x.
    # Remove additive rational terms and logarithms of rational numbers.
    if isinstance(x, sympy.Add):
        terms = x.as_ordered_terms()
    else:
        terms = [x]

    remaining_terms = []
    for term in terms:
        # Rule 1: Remove additive rational terms.
        if term.is_rational:
            continue

        # Rule 2: Remove additive terms that are logarithms of rational numbers.
        # A term 't' is a logarithm of a rational number if exp(t) is rational.
        is_log_of_rational = False
        try:
            val = sympy.exp(term)
            if val.is_rational:
                is_log_of_rational = True
        except (TypeError, ValueError):
            pass  # Cannot determine, assume it's not.
        
        if is_log_of_rational:
            continue
        
        # Keep the term if it's not removed.
        remaining_terms.append(term)
    
    final_answer_expr = sum(remaining_terms)

    # Step 6: Format the final answer in LaTeX.
    latex_answer = sympy.latex(final_answer_expr, ln_notation=True)
    
    # Replace the default log with ln for clarity as per standard notation.
    # The default sympy.latex produces \log{...} for natural logarithm.
    # The ln_notation=True flag helps make it \ln{...}.
    
    print(f"<<<{latex_answer}>>>")

solve_queueing_problem()