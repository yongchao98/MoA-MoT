import sympy

def solve_queueing_problem():
    """
    Solves the queueing theory problem symbolically using SymPy.
    1. Defines the integral for the mean response time x.
    2. Computes the integral.
    3. Identifies and removes terms from x based on the problem's rules.
    4. Formats the final remaining term in LaTeX.
    """
    # Define the symbolic variable
    y = sympy.Symbol('y')

    # Define the arrival rate
    lambda_rate = sympy.Rational(3, 2)

    # Define the job size distribution functions for U(0,1)
    # PDF f_S(y) = 1 for 0 <= y <= 1
    # CDF F_S(y) = y for 0 <= y <= 1
    F_S_y = y

    # Calculate rho_y, the utilization from jobs smaller than y
    # rho_y = lambda * integral(s*f_S(s) ds) from 0 to y
    rho_y = lambda_rate * sympy.integrate(s, (s, 0, y))

    # The integrand for the mean response time x
    integrand = (1 - F_S_y) / (1 - rho_y)

    # Calculate x by integrating from 0 to 1
    # For y > 1, F_S(y) = 1, so the integrand is 0.
    x = sympy.integrate(integrand, (y, 0, 1))

    print("The optimal mean response time is x:")
    # Using sympy.pretty_print for a more readable console output of the expression
    sympy.pretty_print(x, use_unicode=True)
    print("-" * 30)

    # The result x is an 'Add' object, its terms can be accessed via .args
    terms = x.args if isinstance(x, sympy.Add) else [x]
    
    print("Decomposing x into its additive terms:")
    for i, term in enumerate(terms):
        print(f"Term {i+1}:")
        sympy.pretty_print(term, use_unicode=True)
    print("-" * 30)


    def is_rational_log_term(term):
        """
        Checks if a term is a rational multiple of a logarithm of a rational number.
        This corresponds to the interpretation of "additive terms which are logarithms of rational numbers".
        e.g., -4/3 * log(2)
        """
        # A term of the form c, where c is rational, has no log component.
        if term.is_rational:
            return False

        # We are looking for terms of the form: rational_coeff * log(rational_arg)
        # Separate the term into its factors
        factors = term.as_ordered_factors()
        
        log_part = None
        other_coeffs = []
        
        for factor in factors:
            if isinstance(factor, sympy.log):
                log_part = factor
            else:
                other_coeffs.append(factor)

        # Check if there is a log component
        if log_part is None:
            return False

        # Check if the argument of the log is rational
        log_arg = log_part.args[0]
        if not log_arg.is_rational:
            return False
            
        # Check if the combined coefficient is rational
        coeff = sympy.prod(other_coeffs)
        if not coeff.is_rational:
            return False

        return True

    # Identify terms to be removed
    rational_terms_to_remove = [t for t in terms if t.is_rational]
    log_terms_to_remove = [t for t in terms if is_rational_log_term(t)]
    
    terms_to_remove = set(rational_terms_to_remove + log_terms_to_remove)
    
    # Keep the remaining terms
    remaining_terms = [t for t in terms if t not in terms_to_remove]
    
    print("Filtering the terms based on the rules:")
    if not rational_terms_to_remove:
        print("1. No additive rational terms to remove.")
    else:
        print(f"1. Removing rational term(s): {rational_terms_to_remove}")

    if not log_terms_to_remove:
        print("2. No additive 'logarithms of rational numbers' terms to remove.")
    else:
        print(f"2. Removing 'logarithms of rational numbers' term(s):")
        for term in log_terms_to_remove:
             sympy.pretty_print(term, use_unicode=True)

    print("-" * 30)
    
    # There should be a single remaining term as per the problem's phrasing
    if len(remaining_terms) == 1:
        final_term = remaining_terms[0]
        print("The remaining term is:")
        sympy.pretty_print(final_term, use_unicode=True)
        
        # Format the final answer in LaTeX
        # SymPy's default latex output respects the desired order.
        # e.g., (2/3)*sqrt(3)*log(2+sqrt(3)) -> \frac{2 \sqrt{3}}{3} \log{\left(\sqrt{3} + 2 \right)}
        latex_answer = sympy.latex(final_term)
        print("\nFinal Answer in LaTeX format:")
        print(latex_answer)
        
        # The final answer to be submitted
        final_answer_formatted = f"<<<{latex_answer}>>>"
        print("\n" + final_answer_formatted)
        
    else:
        print("Error: Expected one remaining term, but found", len(remaining_terms))
        if remaining_terms:
             print("Remaining terms are:", remaining_terms)
        else:
             print("No terms remaining.")


if __name__ == '__main__':
    solve_queueing_problem()
