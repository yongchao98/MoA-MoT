import sympy

def solve_queueing_problem():
    """
    Solves the queuing theory problem to find the specified remaining term of the
    optimal mean response time x.
    """
    # Step 1: Define symbolic variable and parameters
    s = sympy.Symbol('s')
    lam = sympy.S(3)/2  # Arrival rate
    # Job size S is Uniform(0,1), so its PDF f_S(s) is 1 on [0,1]

    # Step 2: Calculate terms for the E[T(s)] formula
    # rho_s is the traffic intensity from jobs of size <= s
    # rho_s = lam * integral(y * f_S(y) dy) from 0 to s
    # Since f_S(y)=1, integral(y dy) = s**2/2
    rho_s = lam * (s**2 / 2)

    # The busy time contribution from smaller jobs depends on the second moment
    # It is lam * integral(y**2 * f_S(y) dy) from 0 to s
    # Since f_S(y)=1, integral(y**2 dy) = s**3/3
    lam_E_S2_s = lam * (s**3 / 3)

    # Step 3: Define E[T(s)], the mean response time for a job of size s
    # The formula is (lam*E[S_s^2]/2 + s) / (1 - rho_s)
    E_T_s = (lam_E_S2_s / 2 + s) / (1 - rho_s)
    E_T_s_simplified = sympy.simplify(E_T_s)
    
    # Step 4: Integrate E[T(s)] over the job size distribution [0,1] to find x
    # x = integral(E_T_s * f_S(s) ds) from 0 to 1. Since f_S(s)=1, this is:
    x_val = sympy.integrate(E_T_s_simplified, (s, 0, 1))
    
    print("The optimal mean response time, x, is given by the integral:")
    print(f"x = Integral({E_T_s_simplified}, (s, 0, 1))")
    
    # The final equation with each number is printed here
    equation_str = "x = "
    terms_for_print = x_val.as_ordered_terms()
    # Ensure a prettier print order (log term first)
    if terms_for_print[0].is_rational:
         terms_for_print.reverse()
    equation_str += " + ".join(str(t) for t in terms_for_print).replace("+ -", "- ")
    print(f"\nThe value of the integral gives the equation: {equation_str}")


    # Step 5: Analyze and filter the terms of x based on the given rules
    remaining_terms = []
    
    print("\nFiltering the terms of x:")
    for term in x_val.as_ordered_terms():
        # Rule 1: Remove additive rational terms
        if term.is_rational:
            print(f"- The term '{term}' is a rational number. It will be removed.")
            continue
            
        # Rule 2: Remove additive terms which are logarithms of rational numbers
        # A term 't' is a logarithm of a rational number if exp(t) is rational.
        is_log_of_rational = False
        try:
            if sympy.exp(term).is_rational:
                is_log_of_rational = True
        except (TypeError, AttributeError):
            # This handles more complex terms, but is not needed here.
            pass

        if is_log_of_rational:
            print(f"- The term '{term}' is a logarithm of a rational number ({sympy.exp(term)}). It will be removed.")
            continue
        
        # This term is not a logarithm of a rational number. Explain why.
        # This check is more specific for terms of the form r*log(p).
        # We need to check if p^r is rational.
        # e.g., for (16/9)*log(2), check if 2^(16/9) is rational.
        if isinstance(term, sympy.Mul) and any(isinstance(arg, sympy.log) for arg in term.args):
             power = sympy.exp(term)
             print(f"- The term '{term}' equals ln({power}). Since {power} is not rational, this term is not removed.")

        remaining_terms.append(term)
        
    # The final term is the sum of the terms that were not removed
    final_term = sum(remaining_terms)

    # Step 6: Format the answer in LaTeX
    # Manual construction for the desired format "rational multiplicand \ln number"
    if isinstance(final_term, sympy.Mul) and any(isinstance(arg, sympy.log) for arg in final_term.args):
        coeff = final_term.as_coefficient(next(arg for arg in final_term.args if isinstance(arg, sympy.log)))
        log_arg = sympy.exp(final_term/coeff).evalf()
        latex_string = sympy.latex(coeff) + r" \ln " + sympy.latex(sympy.Integer(log_arg))
    else:
        latex_string = sympy.latex(final_term)

    print(f"\nThe remaining term is: {final_term}")
    
    print("\nFinal Answer in LaTeX format:")
    print(latex_string)
    
    # Do not remove the final line.
    print(f"\n<<<{latex_string}>>>")


solve_queueing_problem()